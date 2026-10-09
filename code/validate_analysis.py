"""Independent checks of inputs, arithmetic and delivered outputs; run at package root."""
from pathlib import Path
from itertools import combinations
from io import StringIO
import json
import os
import subprocess
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
rev = root / 'results'
def read(p):
    return pd.read_csv(root / p, sep='\t')
def close(a, b, atol=1e-10, equal_nan=False):
    # Arithmetic equality: no implicit relative tolerance that grows with size.
    assert np.allclose(a, b, rtol=0, atol=atol, equal_nan=equal_nan)
def unique_keys(frame, keys):
    assert not frame[keys].isna().any().any(), ('Missing key', keys)
    assert not frame.duplicated(keys).any(), ('Duplicate key', keys)
def exact_merge(left, right, keys, suffixes=('', '_recalc')):
    # An inner join can hide missing predictions, groups, states or summaries.
    unique_keys(left, keys)
    unique_keys(right, keys)
    joined = left.merge(right, on=keys, how='outer', suffixes=suffixes,
                        validate='one_to_one', indicator=True)
    assert joined['_merge'].eq('both').all(), joined.loc[joined['_merge'] != 'both', keys + ['_merge']]
    return joined.drop(columns='_merge')
def counts_are_valid(values):
    v = np.asarray(values, dtype=float)
    assert np.isfinite(v).all() and (v >= 0).all()
    assert np.equal(v, np.floor(v)).all(), 'Counts must be nonnegative integers'
def bh(values):
    p = np.asarray(values, dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(p)
    out[order] = np.minimum(adjusted, 1)
    return out

checks = {}
manifest = read('data/sample_manifest.tsv')
unique_keys(manifest, ['dataset_id', 'sample_id', 'coarse_cell_type'])
counts_are_valid(manifest.n_cells)
metadata = ['biological_unit_id', 'species', 'age_value', 'age_unit', 'muscle',
            'assay', 'baseline_or_post', 'intervention', 'time_value', 'time_unit']
assert manifest.groupby(['dataset_id', 'sample_id'])[metadata].nunique(dropna=False).le(1).all().all()
assert manifest.baseline_or_post.isin(['baseline', 'post']).all()
checks['manifest_keys_metadata_and_nonnegative_integer_counts']=True

# Independently reconstruct each saved primary task from the manifest captures,
# including its eligible baseline compositions, rather than trust task counts.
tasks = read('data/future_tasks.tsv.gz')
primary = tasks[tasks.primary_task.eq(True)].copy()
unique_keys(primary, ['dataset_id', 'sample_id', 'coarse_cell_type'])
assert len(primary) == 480 and primary.coarse_cell_type.nunique() == 12
assert primary.groupby(['dataset_id', 'sample_id']).ngroups == 40
assert primary.groupby(['dataset_id', 'biological_unit_id']).ngroups == 40
counts_are_valid(primary.n_cells)
counts_are_valid(primary.total_cells)
state_set = set(primary.coarse_cell_type)
base_group = ['biological_unit_id', 'species', 'age_value', 'age_unit', 'muscle',
              'assay', 'baseline_or_post', 'intervention', 'time_value', 'time_unit']
for (dataset, sample), g in primary.groupby(['dataset_id', 'sample_id']):
    assert set(g.coarse_cell_type) == state_set and len(g) == 12
    assert g[['biological_unit_id', 'species', 'age_value', 'age_unit']].nunique(dropna=False).eq(1).all()
    captures = sample.split('+')
    assert len(captures) == len(set(captures))
    rows = manifest[(manifest.dataset_id == dataset) & manifest.sample_id.isin(captures)]
    assert set(rows.sample_id) == set(captures), (dataset, sample)
    assert rows.biological_unit_id.eq(g.biological_unit_id.iloc[0]).all()
    assert rows.baseline_or_post.eq('post').all()
    assert set(rows.coarse_cell_type).issubset(state_set)
    cell_counts = rows.groupby('coarse_cell_type').n_cells.sum().reindex(g.coarse_cell_type, fill_value=0).to_numpy()
    close(g.n_cells, cell_counts, atol=0)
    close(g.total_cells, np.repeat(cell_counts.sum(), len(g)), atol=0)
    close(g.proportion, cell_counts / cell_counts.sum(), atol=1e-12)
    baseline = manifest[(manifest.dataset_id == dataset) & manifest.baseline_or_post.eq('baseline')].copy()
    for field in ['species', 'age_value', 'age_unit']:
        value = g[field].iloc[0]
        baseline = baseline[baseline[field].isna() if pd.isna(value) else baseline[field].eq(value)]
    if dataset == 'human_acute_endurance':
        baseline = baseline[baseline.biological_unit_id.eq(g.biological_unit_id.iloc[0])]
    assert len(baseline), (dataset, sample, 'Missing eligible baseline')
    b = baseline.groupby(base_group + ['coarse_cell_type'], dropna=False).n_cells.sum().unstack('coarse_cell_type', fill_value=0)
    b = b.reindex(columns=g.coarse_cell_type, fill_value=0)
    assert (b.sum(axis=1) > 0).all()
    bp = b.div(b.sum(axis=1), axis=0)
    close(g.baseline_proportion, bp.mean(axis=0), atol=1e-12)
    close(g.baseline_cells, b.sum(axis=0), atol=0)
    close(g.baseline_sd, bp.std(axis=0, ddof=1), atol=1e-12, equal_nan=True)
    close(g.n_baseline_biological_units, np.repeat(len(b), len(g)), atol=0)
checks['all_40_primary_task_counts_and_baselines_reconstructed_from_manifest']=True

# Aggregated annotation counts retain the fixed-threshold input required here;
# individual annotation scores are not supplied or re-thresholded in this package.
a = read('data/annotation_filter_counts.tsv')
unique_keys(a, ['dataset_id', 'sample_id', 'coarse_cell_type'])
counts_are_valid(a.n_before); counts_are_valid(a.n_after)
assert (a.n_after <= a.n_before).all()
assert a.loc[a.coarse_cell_type.eq('other_ambiguous'), 'n_after'].eq(0).all()
rules = {'regeneration_lifespan': 'transfer_score>=0.8',
         'disuse_mechanotherapy': 'marker_margin>=0.3',
         'nmj_growth_atrophy': 'marker_margin>=0.3'}
assert set(a.dataset_id) == set(rules)
assert a.rule.eq(a.dataset_id.map(rules)).all()
declared = manifest[manifest.dataset_id.isin(rules) & manifest.n_cells.gt(0)]
aligned = exact_merge(declared[['dataset_id','sample_id','coarse_cell_type','n_cells']],
                      a[['dataset_id','sample_id','coarse_cell_type','n_before']],
                      ['dataset_id','sample_id','coarse_cell_type'])
close(aligned.n_cells, aligned.n_before, atol=0)
annotation_audit = a.groupby('dataset_id')[['n_before','n_after']].sum().reset_index()
aligned = exact_merge(read('results/annotation_filter_audit.tsv'), annotation_audit, ['dataset_id'])
close(aligned.n_before, aligned.n_before_recalc, atol=0)
close(aligned.n_after, aligned.n_after_recalc, atol=0)
checks['aggregated_annotation_counts_match_manifest_and_filter_totals']=True

# Selective accuracy is recomputed together with coverage, including the N2
# high-score ambiguous predictions; a high score is not calibrated correctness.
transfer = read('data/transfer_predictions.tsv.gz')
assert len(transfer) == 29539 and transfer.held_out_animal.nunique() == 6
assert transfer[['held_out_animal', 'prediction_score', 'predicted', 'truth']].notna().all().all()
assert np.isfinite(transfer.prediction_score).all() and transfer.prediction_score.between(0, 1).all()
selective = []
for animal, g in transfer.groupby('held_out_animal'):
    for threshold in [0, .5, .8]:
        retained = g[g.prediction_score.ge(threshold) & g.predicted.ne('other_ambiguous')]
        selective.append(dict(animal=animal, threshold=threshold, n_total=len(g), n_retained=len(retained),
                              coverage=len(retained)/len(g), accuracy=retained.predicted.eq(retained.truth).mean()))
s = exact_merge(read('results/annotation_selective_validation.tsv'), pd.DataFrame(selective), ['animal', 'threshold'])
for col in ['n_total', 'n_retained', 'coverage', 'accuracy']:
    close(s[col], s[col + '_recalc'], atol=1e-9)
n2 = transfer[transfer.held_out_animal.eq('N2')]
ambiguous_n2 = n2[n2.predicted.eq('other_ambiguous')]
assert len(n2) == 4480 and len(ambiguous_n2) > 0
assert ambiguous_n2.prediction_score.median() > .95
assert round(float(ambiguous_n2.prediction_score.median()), 3) == .957
checks['annotation_selective_accuracy_and_coverage_recalculated']=True
checks['N2_other_ambiguous_score_diagnostic'] = {
    'n_animal': len(n2), 'n_other_ambiguous': len(ambiguous_n2),
    'median_other_ambiguous_score': float(ambiguous_n2.prediction_score.median()),
    'fraction_other_ambiguous_score_at_least_0_8': float(ambiguous_n2.prediction_score.ge(.8).mean())}

p = read('results/benchmark/benchmark_predictions.tsv.gz')
m = read('results/benchmark/benchmark_metrics.tsv')
unique_keys(m, ['split', 'model'])
assert len(m) == 400 and m.model.nunique() == 20 and m.split.nunique() == 20
unique_keys(p, ['split', 'model', 'composition_id', 'coarse_cell_type'])
assert set(p.model) == set(read('model_specifications.tsv').model) == set(m.model)
primary['composition_id'] = primary.dataset_id+'||'+primary.sample_id
split_masks = {}
for field, prefix in [('dataset_id','leave_dataset__'),('intervention','leave_intervention__'),('species','leave_species__')]:
    for value in primary[field].unique():
        split_masks[prefix+str(value).replace(' ','_')] = primary[field].eq(value)
split_masks['time_extrapolation__gt3_days'] = primary.time_days.gt(3)
age_groups = pd.Series(np.where(primary.age_months.isna(),'unknown',np.where(primary.age_months.le(12),'younger','older')),index=primary.index)
for value in age_groups.unique():
    split_masks['leave_age_bin__'+value] = age_groups.eq(value)
assert len(split_masks) == 20 and set(p.split) == set(split_masks) == set(m.split)
for (split, model), prediction_rows in p.groupby(['split','model']):
    required = primary.loc[split_masks[split],['composition_id','coarse_cell_type','proportion','baseline_proportion','n_cells']]
    aligned = exact_merge(prediction_rows[['composition_id','coarse_cell_type','proportion','baseline_proportion','n_cells']],
                          required,['composition_id','coarse_cell_type'])
    for col in ['proportion','baseline_proportion','n_cells']:
        close(aligned[col],aligned[col+'_recalc'],atol=1e-12)
checks['every_model_split_has_exactly_the_required_input_targets'] = True
assert np.isfinite(p.prediction).all() and p.prediction.between(0, 1).all()
assert p.groupby(['split', 'model', 'composition_id']).coarse_cell_type.nunique().eq(12).all()
close(p.groupby(['split', 'model', 'composition_id']).prediction.sum(), 1, atol=1e-10)
calc = p.assign(se=(p.prediction-p.proportion)**2, ae=abs(p.prediction-p.proportion)).groupby(['split', 'model']).agg(se=('se','mean'), mae=('ae','mean')).reset_index()
check = exact_merge(m, calc, ['split','model'])
close(check.rmse, np.sqrt(check.se), atol=1e-12)
close(check.mae, check.mae_recalc, atol=1e-12)
checks['400_model_split_errors_independently_recalculated'] = True
js = []
for key, g in p.groupby(['split', 'model', 'composition_id'], sort=False):
    a = np.maximum(g.proportion.to_numpy(), 1e-12); a /= a.sum()
    b = np.maximum(g.prediction.to_numpy(), 1e-12); b /= b.sum(); mid = (a+b)/2
    js.append((*key, .5*np.sum(a*np.log2(a/mid))+.5*np.sum(b*np.log2(b/mid))))
j = pd.DataFrame(js, columns=['split','model','composition_id','js']).groupby(['split','model']).js.mean().reset_index()
z = exact_merge(m, j, ['split','model'])
close(z.mean_js_divergence, z.js, atol=1e-12)
checks['JS_independently_recalculated'] = True
conv = read('results/benchmark/model_convergence_audit.tsv')
unique_keys(conv, ['split'])
assert set(conv.split) == set(m.split) and len(conv) == 20 and conv.convergence.eq(0).all()
checks['all_20_neural_network_fits_converged'] = True
for sub in ['benchmark', 'controls_exclude_ambiguous']:
    f = read(f'results/{sub}/negative_control_split_results.tsv')
    sm = read(f'results/{sub}/negative_control_summary.tsv')
    unique_keys(f,['held_out_dataset','control','replicate'])
    assert set(f.held_out_dataset) == set(primary.dataset_id)
    for (dataset, control), group in f.groupby(['held_out_dataset','control']):
        ids = set(range(1,26)) if control in {'baseline_vector_permutation','intervention_label_permutation','timepoint_permutation'} else {0}
        assert set(group.replicate) == ids
    assert np.isfinite(f.rmse).all() and f.rmse.ge(0).all()
    means = f.groupby(['held_out_dataset','control']).rmse.mean().reset_index()
    assert means.groupby('control').held_out_dataset.nunique().eq(5).all()
    average = means.groupby('control').rmse.mean().rename('recalculated_mean').reset_index()
    aligned = exact_merge(sm, average, ['control'])
    close(aligned.mean_rmse, aligned.recalculated_mean, atol=1e-12)
    paired = read(f'results/{sub}/negative_control_paired_comparisons.tsv')
    unique_keys(paired, ['control'])
    assert set(paired.control) == set(means.control)-{'observed_baseline_context'}
    assert not any('signflip' in c or c.endswith('_p') for c in paired.columns)
    observed = means[means.control.eq('observed_baseline_context')][['held_out_dataset','rmse']]
    for _, g in paired.iterrows():
        aligned = exact_merge(means[means.control.eq(g.control)][['held_out_dataset','rmse']], observed,
                              ['held_out_dataset'], suffixes=('_control','_observed'))
        difference = aligned.rmse_control-aligned.rmse_observed
        close(g.mean_rmse_advantage_of_observed, difference.mean(), atol=1e-12)
        assert g.observed_better_count == (difference > 0).sum()
    checks[sub+'_control_arithmetic_and_no_invalid_fold_P_values'] = True

t = read('results/tasks_exclude_ambiguous.tsv.gz')
unique_keys(t, ['dataset_id','sample_id','coarse_cell_type'])
assert len(t) == 440 and t.coarse_cell_type.nunique() == 11
assert not t.coarse_cell_type.eq('other_ambiguous').any()
counts_are_valid(t.n_cells); counts_are_valid(t.total_cells)
close(t.n_cells/t.total_cells, t.proportion, atol=1e-12)
close(t.groupby(['dataset_id','sample_id']).n_cells.sum(), t.groupby(['dataset_id','sample_id']).total_cells.first(), atol=0)
close(t.groupby(['dataset_id','sample_id']).baseline_proportion.sum(), 1, atol=1e-12)
assert (t.capture_total_cells >= t.total_cells).all()
checks['11_category_counts_denominators_and_baselines_consistent'] = True

# All 15 participant assignments and BH adjustment, independent of the R test.
h = read('data/human_paired_changes.tsv')
unique_keys(h, ['subject','coarse_cell_type'])
assert h.groupby('subject').Group.nunique().eq(1).all()
subjects = h[['subject','Group']].drop_duplicates().sort_values('subject')
assert len(subjects) == 6 and subjects.Group.eq('Exercise').sum() == 4 and subjects.Group.eq('Control').sum() == 2
assert h.groupby('coarse_cell_type').subject.nunique().eq(6).all()
assert h[['Pre','Post']].ge(0).all().all() and h[['Pre','Post']].le(1).all().all()
close(h.Post-h.Pre, h.post_minus_pre, atol=1e-12)
close(h.groupby('subject')[['Pre','Post']].sum(), 1, atol=1e-12)
assignments = list(combinations(subjects.subject, 4)); assert len(assignments) == 15
human = []
for state, g in h.groupby('coarse_cell_type'):
    ex = g.loc[g.Group.eq('Exercise'),'post_minus_pre'].mean()
    control = g.loc[g.Group.eq('Control'),'post_minus_pre'].mean()
    observed = ex-control
    null = [g.loc[g.subject.isin(ids),'post_minus_pre'].mean()-g.loc[~g.subject.isin(ids),'post_minus_pre'].mean() for ids in assignments]
    human.append(dict(coarse_cell_type=state, exercise_mean_paired_change=ex, control_mean_paired_change=control,
                      difference_in_paired_change=observed, exact_two_sided_p=np.mean(np.abs(null)>=abs(observed)-1e-15)))
human = pd.DataFrame(human); human['fdr_bh'] = bh(human.exact_two_sided_p)
aligned = exact_merge(read('results/human_exact_tests.tsv'), human, ['coarse_cell_type'])
for col in ['exercise_mean_paired_change','control_mean_paired_change','difference_in_paired_change','exact_two_sided_p','fdr_bh']:
    close(aligned[col], aligned[col+'_recalc'], atol=1e-12)
assert aligned.n_exact_assignments.eq(15).all()
assert (aligned.exact_two_sided_p >= 1/15-1e-12).all()
checks['human_15_assignments_and_BH_independently_recalculated'] = True

# Reconstruct both denervation representations and all 20 animal assignments.
den = read('data/denervation_published_counts.tsv')
unique_keys(den, ['sample_id','published_cell_type'])
counts_are_valid(den.n_cells)
assert len(den) == 78 and den.sample_id.nunique() == 6 and den.published_cell_type.nunique() == 13
assert den.n_cells.sum() == 29539
assert den.groupby('sample_id').condition.nunique().eq(1).all()
assert set(den.loc[den.condition.eq('control'),'sample_id']) == {'GSM5571584','GSM5571585','GSM5571586'}
assert set(den.loc[den.condition.eq('denervation'),'sample_id']) == {'GSM5571587','GSM5571588','GSM5571589'}
label_map = {'TypeI':'myonuclei','TypeIIA':'myonuclei','TypeIIx':'myonuclei','TypeIIB':'myonuclei',
             'MTJ':'tendon_mtj','NMJ':'nmj_myonuclei','Satellite':'muscle_stem_myogenic',
             'FAPs':'fap_mesenchymal','Pericyte':'pericyte','Adipocyte':'adipocyte',
             'Endothelial':'endothelial','Macrophage':'immune','Unclassified':'other_ambiguous'}
assert set(den.published_cell_type) == set(label_map)
den['coarse_cell_type'] = den.published_cell_type.map(label_map)
coarse = den.groupby(['sample_id','condition','coarse_cell_type']).n_cells.sum().reset_index()
aligned = exact_merge(coarse[['sample_id','coarse_cell_type','n_cells']],
    manifest[manifest.dataset_id.eq('denervation')][['sample_id','coarse_cell_type','n_cells']],
    ['sample_id','coarse_cell_type'])
close(aligned.n_cells, aligned.n_cells_recalc, atol=0)
assert len(coarse) == 60 and coarse.coarse_cell_type.nunique() == 10
animal_ids = sorted(den.sample_id.unique())
combinations_20 = list(combinations(animal_ids,3)); assert len(combinations_20) == 20
output_composition = read('results/denervation_animal_composition.tsv')
output_tests = read('results/denervation_composition_exact_tests.tsv')
unique_keys(output_composition,['representation','sample_id','state'])
unique_keys(output_tests,['representation','state'])
assert len(output_tests) == 23 and len(output_composition) == 138
assert set(output_tests.representation) == {'published_13_labels','observed_coarse_10_categories'}
for representation, data, state_column, n_states in [
    ('published_13_labels',den,'published_cell_type',13),
    ('observed_coarse_10_categories',coarse,'coarse_cell_type',10)]:
    data = data.rename(columns={state_column:'state'}).copy()
    data['total_cells'] = data.groupby('sample_id').n_cells.transform('sum')
    data['proportion'] = data.n_cells/data.total_cells
    aligned = exact_merge(output_composition[output_composition.representation.eq(representation)],
        data[['sample_id','state','n_cells','total_cells','proportion']],['sample_id','state'])
    close(aligned.n_cells,aligned.n_cells_recalc,atol=0)
    close(aligned.total_cells,aligned.total_cells_recalc,atol=0)
    close(aligned.proportion,aligned.proportion_recalc,atol=1e-12)
    tests = []
    for state, g in data.groupby('state'):
        assert len(g) == 6 and set(g.sample_id) == set(animal_ids)
        control = g.loc[g.condition.eq('control'),'proportion'].mean()
        changed = g.loc[g.condition.eq('denervation'),'proportion'].mean()
        difference = changed-control
        null = [g.loc[g.sample_id.isin(ids),'proportion'].mean()-g.loc[~g.sample_id.isin(ids),'proportion'].mean() for ids in combinations_20]
        tests.append(dict(state=state,control_mean=control,denervation_mean=changed,
                          denervation_minus_control=difference,
                          exact_two_sided_p=np.mean(np.abs(null)>=abs(difference)-1e-15)))
    tests = pd.DataFrame(tests); tests['fdr_bh'] = bh(tests.exact_two_sided_p)
    assert len(tests) == n_states
    aligned = exact_merge(output_tests[output_tests.representation.eq(representation)],tests,['state'])
    for col in ['control_mean','denervation_mean','denervation_minus_control','exact_two_sided_p','fdr_bh']:
        close(aligned[col],aligned[col+'_recalc'],atol=1e-12)
    assert aligned.n_assignments.eq(20).all()
    close(aligned.minimum_attainable_two_sided_p,.1,atol=1e-12)
    assert (aligned.exact_two_sided_p >= .1-1e-12).all() and aligned.fdr_bh.ge(.05).all()
    close(aligned.exact_two_sided_p.min(),.1 if n_states == 13 else .2,atol=1e-12)
checks['denervation_13_fine_and_10_coarse_tests_recomputed_from_counts'] = True
de = read('data/denervation_gene_statistics.tsv.xz')
unique_keys(de,['cell_type','gene'])
assert de.cell_type.nunique() == 7 and de.padj.dropna().between(0,1).all()
assert np.isfinite(de.log2FoldChange).all()
de_counts = de.groupby('cell_type').apply(lambda g: pd.Series({'n_gene_rows':len(g),
    'n_nonmissing_padj':g.padj.notna().sum(),'n_fdr_005':g.padj.lt(.05).sum(),
    'n_fdr_005_abs_lfc_1':(g.padj.lt(.05)&g.log2FoldChange.abs().ge(1)).sum()}),include_groups=False).reset_index()
aligned = exact_merge(read('results/denervation_pseudobulk_summary.tsv'),de_counts,['cell_type'])
for col in ['n_gene_rows','n_nonmissing_padj','n_fdr_005','n_fdr_005_abs_lfc_1']:
    close(aligned[col],aligned[col+'_recalc'],atol=0)
assert de_counts.set_index('cell_type').n_fdr_005.to_dict() == {
    'TypeIIB':3169,'FAPs':567,'TypeIIA':420,'TypeIIx':279,'Macrophage':126,'Endothelial':124,'Unclassified':811}
checks['denervation_DE_counts_recomputed_from_supplied_gene_statistics'] = True

spatial = read('data/spatial_gene_statistics.tsv.xz')
unique_keys(spatial, ['gene'])
g = spatial[spatial.padj < .05]
assert len(g) == 907
same = np.sign(g.log2FoldChange) == np.sign(g.log2FoldChange_excluding_PH078)
assert same.sum() == 907 and (same & (g.padj_excluding_PH078 < .05)).sum() == 103
close(g.exact_signflip_p.min(), 2/2**5, atol=1e-12)
checks['spatial_907_directions_103_significant_and_P_resolution'] = True
for source, out, replicates in [('baseline_resampling_unit_draws.tsv.gz','baseline_resampling_group_summary.tsv',2000),
                                ('rarefaction_unit_draws.tsv.gz','rarefaction_summary.tsv',200)]:
    d = read('results/'+source)
    unique_keys(d, ['id','replicate'])
    assert d.groupby('id').replicate.nunique().eq(replicates).all()
    assert set(d.replicate) == set(range(replicates))
    assert np.isfinite(d.js).all() and d.js.between(0,1).all()
    group = d.groupby(['dataset_id','intervention','replicate']).js.median()
    q = group.groupby(['dataset_id','intervention']).agg(median='median', q025=lambda x:x.quantile(.025), q975=lambda x:x.quantile(.975)).reset_index()
    v = exact_merge(read('results/'+out), q, ['dataset_id','intervention'])
    for col in ['median','q025','q975']:
        close(v[col], v[col+'_recalc'], atol=1e-9)
    checks[source+'_summary_recalculated'] = True

# Audit the count input and model design independently of the fitting script.
# R is already a required runtime; this read-only check requires no new package.
r_audit = r'''
x <- readRDS("data/lineage_counts.rds"); m <- x$sample_metadata; counts <- x$collapsed_counts
stopifnot(identical(colnames(counts),as.character(m$sample)),!anyDuplicated(m$sample),
 !anyDuplicated(rownames(counts)),all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
cat("contrast\tn_samples\tn_parameters\tdesign_rank\tn_genes\n")
def <- list(c("DEV_POS","DEV_NEG","developmental_GFP_positive_vs_negative"),
 c("MOV_POS","MOV_NEG","overload_newly_fused_GFP_positive_vs_preexisting_GFP_negative"),
 c("MYMK_KO","SHAM","overload_fusion_blockade_vs_sham_batch_adjusted"))
for(k in seq_along(def)) {
 d <- def[[k]]; md <- droplevels(m[m$condition %in% d[1:2],,drop=FALSE]); md$condition<-factor(md$condition)
 md$preparation_pair<-factor(md$preparation_pair);md$batch<-factor(md$batch)
 if(k<3L) stopifnot(all(table(md$preparation_pair,md$condition)==1L),nlevels(md$preparation_pair)==3L)
 X <- if(k<3L) model.matrix(~preparation_pair+condition,md) else model.matrix(~batch+condition,md)
 nc <- sum(rowSums(counts[,md$sample,drop=FALSE]>=10)>=2)
 cat(paste(d[3],nrow(md),ncol(X),qr(X)$rank,nc,sep="\t"),"\n",sep="")
}
'''
r_output = subprocess.run([os.environ.get('RSCRIPT','Rscript'),'-e',r_audit],cwd=root,check=True,capture_output=True,text=True).stdout
count_audit = pd.read_csv(StringIO(r_output),sep='\t')
lineage_audit = exact_merge(read('results/lineage_design_audit.tsv'), count_audit, ['contrast'])
assert len(lineage_audit) == 3
for col in ['n_samples','n_parameters','design_rank','n_genes']:
    close(lineage_audit[col], lineage_audit[col+'_recalc'], atol=0)
assert lineage_audit.design_rank.eq(lineage_audit.n_parameters).all()
assert lineage_audit.final_failed_coefficients.eq(0).all()
lineage = read('results/lineage_gene_statistics.tsv')
unique_keys(lineage, ['contrast','gene'])
assert set(lineage.contrast) == set(lineage_audit.contrast)
assert lineage.padj.dropna().between(0,1).all()
assert lineage.pvalue.dropna().between(0,1).all()
assert (lineage.lfcSE.dropna() >= 0).all()
summary = lineage.groupby('contrast').apply(lambda g: pd.Series({'n_genes_tested':len(g),
    'n_fdr_005':g.padj.lt(.05).sum(),'n_fdr_005_abs_lfc_1':(g.padj.lt(.05)&g.log2FoldChange.abs().ge(1)).sum()}),include_groups=False).reset_index()
aligned = exact_merge(read('results/gold_lineage_summary.tsv'), summary, ['contrast'])
for col in ['n_genes_tested','n_fdr_005','n_fdr_005_abs_lfc_1']:
    close(aligned[col], aligned[col+'_recalc'], atol=0)
for name, significant, large in [('developmental_GFP_positive_vs_negative',10173,6015),
    ('overload_newly_fused_GFP_positive_vs_preexisting_GFP_negative',441,49),
    ('overload_fusion_blockade_vs_sham_batch_adjusted',1219,753)]:
    row = summary.set_index('contrast').loc[name]
    assert row.n_fdr_005 == significant and row.n_fdr_005_abs_lfc_1 == large
close(summary.set_index('contrast').n_genes_tested.sort_index(), lineage_audit.set_index('contrast').n_genes.sort_index(),atol=0)
checks['lineage_count_input_design_ranks_gene_counts_and_all_coefficient_fits_converged'] = True
checks['lineage_results_refitted_from_counts'] = True

qc = read('data/reference_qc_overlap.tsv')
assert qc.n_publication_retained.sum() == 29539
assert qc.n_our_qc_singlet_pass.sum() < qc.n_publication_retained.sum()
checks['discovery_accounting_and_publication_vs_local_QC_distinction'] = True
assert len(list((root/'figures').glob('Figure*.png'))) == 5
checks['five_figures_generated'] = True
(rev/'analysis_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
