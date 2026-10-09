#!/usr/bin/env python3
"""R1.1, R1.2, R1.4: reproducible sensitivity analyses; run from project root."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; OUT.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(20261008)
def read(p): return pd.read_csv(ROOT/p,sep='\t')
def save(x,n): x.to_csv(OUT/n,sep='\t',index=False,float_format='%.10g')
m=read('data/sample_manifest.tsv'); m=m[m.dataset_id!='fusion_population'].copy()
fields=['dataset_id','biological_unit_id','species','age_value','age_unit','muscle','assay','baseline_or_post','intervention','time_value','time_unit']
m[fields]=m[fields].fillna('unspecified').astype(str)
meta=m.groupby(fields,dropna=False).sample_id.agg(lambda v:'+'.join(sorted(set(v)))).reset_index()
meta['id']=meta.dataset_id+'||'+meta.sample_id
m=m.merge(meta[fields+['id']],on=fields,validate='many_to_one')
counts=m.pivot_table(index='id',columns='coarse_cell_type',values='n_cells',aggfunc='sum',fill_value=0)
# Complete the exact archived 12-state universe, including zero counts where a state is absent from a unit.
tasks=read('data/future_tasks.tsv.gz'); states=sorted(tasks.coarse_cell_type.unique())
counts=counts.reindex(columns=states,fill_value=0); meta=meta.set_index('id').loc[counts.index]
meta['key']=meta.dataset_id+'||'+meta.species+'||'+meta.age_value+'||'+meta.age_unit
paired=meta.dataset_id=='human_acute_endurance'
meta.loc[paired,'key']+='||'+meta.loc[paired,'biological_unit_id']
assert len(states)==12

def js(p,q):
 p=np.asarray(p,float);q=np.asarray(q,float);p=p/p.sum(axis=-1,keepdims=True);q=q/q.sum(axis=-1,keepdims=True);mid=(p+q)/2
 with np.errstate(divide='ignore',invalid='ignore'):
  return .5*np.sum(np.where(p>0,p*np.log2(p/mid),0),axis=-1)+.5*np.sum(np.where(q>0,q*np.log2(q/mid),0),axis=-1)

def evaluate(c,scenario,subset=None,min_total=1):
 md=meta.loc[c.index].copy(); keep=c.sum(axis=1)>=min_total
 if subset is not None: keep &= subset.reindex(c.index,fill_value=False)
 md=md[keep];c=c.loc[md.index];p=c.div(c.sum(axis=1),axis=0);rows=[]
 for idx,z in md[md.baseline_or_post=='post'].iterrows():
  b=md[(md.key==z.key)&(md.baseline_or_post=='baseline')].index
  if not len(b):continue
  a=p.loc[b].to_numpy();q=p.loc[idx].to_numpy();v=js(a.mean(0),q)
  row=z.to_dict();row.update(id=idx,scenario=scenario,n_baselines=len(b),n_future_cells=int(c.loc[idx].sum()),js=v,
    individual_baseline_min=float(js(a,q).min()),individual_baseline_median=float(np.median(js(a,q))),individual_baseline_max=float(js(a,q).max()))
  rows.append(row)
 return pd.DataFrame(rows)
allres=[evaluate(counts,'archived_12_states')]
primary=allres[0]
arch=read('results/future_unit_dynamics.tsv')
check=primary.merge(arch,on=['dataset_id','biological_unit_id'],suffixes=('','_arch'))
assert len(check)==40 and np.max(abs(check.js-check.plasticity_index_js))<1e-9
known=counts.drop(columns='other_ambiguous')
for name,c,sub,minimum in [
 ('exclude_ambiguous',known,None,1),
 ('exclude_ambiguous_min2000',known,None,2000),
 ('mouse_only',counts,meta.species=='Mus musculus',1),
 ('snRNA_only',counts,meta.assay=='snRNA-seq',1),
 ('mouse_snRNA_only',counts,(meta.species=='Mus musculus')&(meta.assay=='snRNA-seq'),1),
 ('exclude_N2',counts,meta.biological_unit_id!='GSM5571588',1),
 ('exclude_ambiguous_and_N2',known,meta.biological_unit_id!='GSM5571588',1)]:
 allres.append(evaluate(c,name,sub,minimum))
# Fixed-threshold per-sample/state counts are supplied as processed inputs.
filtered=m.copy(); audit=[]
annotation_counts=read('data/annotation_filter_counts.tsv')
for d,a in annotation_counts.groupby('dataset_id',sort=False):
 cnt=a.set_index(['sample_id','coarse_cell_type']).n_after
 sel=filtered.dataset_id==d
 filtered.loc[sel,'n_cells']=[cnt.get((sample,state),0) for sample,state in zip(filtered.loc[sel,'sample_id'],filtered.loc[sel,'coarse_cell_type'])]
 audit.append(dict(dataset_id=d,rule=a.rule.iloc[0],n_before=int(a.n_before.sum()),n_after=int(a.n_after.sum())))
f=filtered.pivot_table(index='id',columns='coarse_cell_type',values='n_cells',aggfunc='sum',fill_value=0).reindex(index=counts.index,columns=states,fill_value=0).drop(columns='other_ambiguous')
allres.append(evaluate(f,'higher_annotation_confidence'))
save(pd.DataFrame(audit),'annotation_filter_audit.tsv')
# Baseline resampling retains future compositions and resamples full analysis-unit vectors.
boot=[];basepairs=[];pp=counts.div(counts.sum(axis=1),axis=0)
for key,g in meta.groupby('key'):
 b=g[g.baseline_or_post=='baseline'].index; fut=g[g.baseline_or_post=='post'].index
 if not len(b) or not len(fut):continue
 a=pp.loc[b].to_numpy(); draws=rng.integers(0,len(b),size=(2000,len(b)));means=a[draws].mean(axis=1)
 for i in fut:
  v=js(means,pp.loc[i].to_numpy())
  for rep,val in enumerate(v):boot.append(dict(id=i,dataset_id=meta.loc[i,'dataset_id'],intervention=meta.loc[i,'intervention'],replicate=rep,js=val,n_baselines=len(b)))
 for j in range(len(b)):
  for k in range(j):basepairs.append(dict(key=key,baseline_1=b[j],baseline_2=b[k],js=float(js(a[j],a[k]))))
bt=pd.DataFrame(boot);bs=bt.groupby(['dataset_id','intervention','replicate']).js.median().reset_index()
sumboot=bs.groupby(['dataset_id','intervention']).js.agg(median='median',q025=lambda x:x.quantile(.025),q975=lambda x:x.quantile(.975)).reset_index()
save(bt,'baseline_resampling_unit_draws.tsv.gz');save(bs,'baseline_resampling_group_draws.tsv.gz')
save(sumboot,'baseline_resampling_group_summary.tsv');save(pd.DataFrame(basepairs),'baseline_pairwise_js.tsv')
save(primary[['id','dataset_id','intervention','n_baselines','js','individual_baseline_min','individual_baseline_median','individual_baseline_max']],'individual_baseline_sensitivity.tsv')
# Equal-depth rarefaction at 500 observed assigned cells, without replacement.
rares=[];eligible=known[known.sum(axis=1)>=500].astype(int)
for rep in range(200):
 rc=pd.DataFrame([rng.multivariate_hypergeometric(v,500) for v in eligible.to_numpy()],index=eligible.index,columns=eligible.columns)
 z=evaluate(rc,'rarefaction_500');z['replicate']=rep;rares.append(z)
rr=pd.concat(rares,ignore_index=True);save(rr,'rarefaction_unit_draws.tsv.gz');save(meta.assign(assigned_cells=known.sum(axis=1),eligible_500=known.sum(axis=1)>=500).reset_index(),'rarefaction_eligibility.tsv');rmed=rr.groupby(['dataset_id','intervention','replicate']).js.median().reset_index()
save(rmed.groupby(['dataset_id','intervention']).js.agg(median='median',q025=lambda x:x.quantile(.025),q975=lambda x:x.quantile(.975)).reset_index(),'rarefaction_summary.tsv')
res=pd.concat(allres,ignore_index=True);save(res,'occupancy_sensitivity_units.tsv')
summary=res.groupby(['scenario','dataset_id','intervention']).agg(n_future=('id','nunique'),median_js=('js','median'),q25=('js',lambda x:x.quantile(.25)),q75=('js',lambda x:x.quantile(.75)),minimum_baselines=('n_baselines','min'),maximum_baselines=('n_baselines','max')).reset_index();save(summary,'occupancy_sensitivity_summary.tsv')
# Descriptive subset inventory; no re-estimation of within-study distances. Omit one dataset at a time; no implied replacement for an omitted perturbation.
loo=[]
for omit in sorted(primary.dataset_id.unique()):
 z=primary[primary.dataset_id!=omit].groupby(['dataset_id','intervention']).js.median().sort_values(ascending=False).reset_index();z['excluded_dataset']=omit;z['rank_remaining']=range(1,len(z)+1);loo.append(z)
save(pd.concat(loo),'leave_dataset_out_ordering.tsv')
# Confidence threshold cannot diagnose misclassification with high confidence.
pred=read('data/transfer_predictions.tsv.gz')
confusion=pred.groupby(['held_out_animal','truth','predicted']).size().rename('n_cells').reset_index()
save(confusion,'annotation_confusion.tsv')
ambiguity=pred[pred.predicted=='other_ambiguous'].groupby('held_out_animal').agg(
 n_predicted_ambiguous=('predicted','size'),mean_score=('prediction_score','mean'),
 median_score=('prediction_score','median')).reset_index()
save(ambiguity,'annotation_ambiguous_scores.tsv')
a=[]
for animal,g in pred.groupby('held_out_animal'):
 for threshold in [0,.5,.8]:
  sub=g[(g.prediction_score>=threshold)&(g.predicted!='other_ambiguous')]
  a.append(dict(animal=animal,threshold=threshold,n_total=len(g),n_retained=len(sub),coverage=len(sub)/len(g),accuracy=(sub.predicted==sub.truth).mean()))
save(pd.DataFrame(a),'annotation_selective_validation.tsv')
# Audit age imputation and sample/state accounting directly from saved primary tasks.
t=tasks[tasks.primary_task].copy();save(pd.DataFrame([dict(n_rows=len(t),n_samples=t.groupby(['dataset_id','sample_id']).ngroups,n_biological_units=t.groupby(['dataset_id','biological_unit_id']).ngroups,n_states=t.coarse_cell_type.nunique(),missing_age_rows=t.age_months.isna().sum(),all_task_age_median=t.age_months.median())]),'primary_task_audit.tsv')
print(summary.to_string(index=False));print('\nBASELINE RESAMPLING\n',sumboot.to_string(index=False))
# Export a coherent 11-state task table: renormalize EACH eligible baseline before averaging.
kprop=known.div(known.sum(axis=1),axis=0)
kt=tasks[tasks.primary_task & (tasks.coarse_cell_type!='other_ambiguous')].copy()
kt['capture_total_cells']=kt.total_cells
kt['capture_n_cells']=kt.n_cells
for ix,row in kt.iterrows():
 uid=row.dataset_id+'||'+row.sample_id; state=row.coarse_cell_type
 base=meta[(meta.key==meta.loc[uid,'key'])&(meta.baseline_or_post=='baseline')].index
 kt.at[ix,'n_cells']=known.loc[uid,state]
 kt.at[ix,'total_cells']=known.loc[uid].sum()
 kt.at[ix,'baseline_cells']=known.loc[base,state].sum()
 kt.at[ix,'baseline_sd']=kprop.loc[base,state].std(ddof=1)
 kt.at[ix,'proportion']=kprop.loc[uid,state]
 kt.at[ix,'baseline_proportion']=kprop.loc[base,state].mean()
kt.to_csv(OUT/'tasks_exclude_ambiguous.tsv.gz',sep='\t',index=False,na_rep='NA')
