#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("code/io.R")
source("code/style.R")
theme_set(theme_nature_contract(base_family="sans"))


fig_dir <- "figures"
src_dir <- "results/figure_source_data"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
dataset_labels <- c(fusion_population="Resting human", regeneration_lifespan="Regeneration", nmj_growth_atrophy="NMJ atlas", disuse_mechanotherapy="Disuse", human_acute_endurance="Acute exercise", denervation="Denervation")
label_dataset <- function(x) {
 out <- unname(dataset_labels[x]); out[is.na(out)] <- gsub("_", " ", x[is.na(out)]); out
}
label_intervention <- function(x) {
 labels <- c(notexin_injury="Notexin injury", mechanical_overload="Overload", hindlimb_unloading_reloading="Unloading / reloading", hindlimb_unloading_reloading_mechanotherapy="Reloading + mechanotherapy", circadian_rest_control="Circadian control", acute_endurance_exercise="Acute exercise", sciatic_denervation="Sciatic denervation", denervation="Denervation", basal_control="Basal control")
 out <- unname(labels[x]); out[is.na(out)] <- gsub("_", " ", x[is.na(out)]); out
}
label_split <- function(x) {
 out <- gsub("leave_dataset__", "Study: ", x)
 for(n in names(dataset_labels)) out <- gsub(n, dataset_labels[[n]], out, fixed=TRUE)
 out <- gsub("leave_intervention__", "Intervention: ", out)
 out <- gsub("hindlimb_unloading_reloading_mechanotherapy", "reloading + mech.", out, fixed=TRUE)
 out <- gsub("hindlimb_unloading_reloading", "unloading / reloading", out, fixed=TRUE)
 out <- gsub("acute_endurance_exercise", "acute exercise", out, fixed=TRUE)
 out <- gsub("circadian_rest_control", "circadian control", out, fixed=TRUE)
 out <- gsub("leave_species__", "Species: ", out)
 out <- gsub("Mus_musculus", "mouse", out, fixed=TRUE)
 out <- gsub("Homo_sapiens", "human", out, fixed=TRUE)
 out <- gsub("Rattus_norvegicus", "rat", out, fixed=TRUE)
 out <- gsub("leave_age_bin__", "Age: ", out)
 out <- gsub("time_extrapolation__gt3_days", "Time: >3 days", out, fixed=TRUE)
 gsub("_", " ", out)
}

# Figure 1: scope, thresholds, QC, evidence.
manifest <- fread("data/sample_manifest.tsv")
scope <- manifest[, .(n_cells=sum(n_cells),n_observed_combinations=.N,n_biological_units=uniqueN(biological_unit_id)),by=dataset_id]
scope <- merge(scope,fread("data/registry.tsv"),by="dataset_id",all.x=TRUE)
scope[, dataset_label := label_dataset(dataset_id)]
atomic_write_tsv(scope, file.path(src_dir, "figure1_scope.tsv"))
p1a <- ggplot(scope, aes(reorder(dataset_label, n_cells), n_cells / 1000, fill = evidence_tier)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(Baseline = unname(palette_contract["neutral_mid"]), Bronze = "#BDBDBD", Silver = "#80B1D3", Gold = unname(palette_contract["gold"]))) +
  scale_y_continuous(breaks=c(0,100,200)) +
  labs(x = NULL, y = "Cells/nuclei (thousands)", fill = NULL) +
  guides(fill=guide_legend(nrow=2)) + theme(legend.position="bottom")
p1b <- ggplot(scope, aes(reorder(dataset_label, n_observed_combinations), n_observed_combinations, fill = evidence_tier)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(Baseline = unname(palette_contract["neutral_mid"]), Bronze = "#BDBDBD", Silver = "#80B1D3", Gold = unname(palette_contract["gold"])), guide = "none") +
  scale_y_continuous(breaks=c(0,1000,2000), labels=scales::label_number(scale_cut=scales::cut_short_scale())) +
  labs(x = NULL, y = "Sample–state combinations")
qc <- fread("data/qc.tsv")
qc[, retained_fraction := n_final_pass / n_input]
atomic_write_tsv(qc, file.path(src_dir, "figure1_qc.tsv"))
p1c <- ggplot(qc, aes(dataset, retained_fraction, colour = dataset)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35) + geom_point(position=position_jitter(width=0.12,height=0,seed=20261008),size=0.8,alpha=0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Dark2", guide = "none") + labs(x = NULL, y = "QC retention") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p1d <- ggplot(scope, aes(evidence_tier, n_biological_units, size = n_cells, colour = evidence_tier)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c(Baseline = unname(palette_contract["neutral_mid"]), Bronze = "#8C8C8C", Silver = "#377EB8", Gold = unname(palette_contract["gold"])), guide = "none") +
  scale_y_continuous(expand=expansion(mult=c(0.08,0.18))) +
  scale_size_area(max_size = 7,breaks=c(50000,100000,150000),labels=c("50k","100k","150k")) + labs(x = "Evidence tier", y = "Analysis units", size = "Cells/nuclei") +
  theme(legend.position="bottom")
fig1 <- (p1a | p1b) / (p1c | p1d) + plot_annotation(tag_levels = "a")


# Figure 5: strict generalization, controls, and in-sample variance decomposition.
metrics <- fread("results/benchmark/benchmark_metrics.tsv")
selected_splits <- metrics[grepl("leave_dataset|leave_species|time_extrapolation", split)]
top_models <- metrics[, .(score = mean(mean_js_divergence)), by = model][order(score)][1:8, model]
heat <- selected_splits[model %in% top_models]
heat[, split_label := label_split(split)]
atomic_write_tsv(heat, file.path(src_dir, "figure5_benchmark_heatmap.tsv"))
p2a <- ggplot(heat, aes(model, split_label, fill = mean_js_divergence)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  labs(x = NULL, y = NULL, fill = "Mean JS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
controls <- fread("results/benchmark/negative_control_summary.tsv")
atomic_write_tsv(controls, file.path(src_dir, "figure5_negative_controls.tsv"))
control_pairs <- fread("results/benchmark/negative_control_paired_comparisons.tsv")
atomic_write_tsv(control_pairs, file.path(src_dir, "figure5_negative_control_paired.tsv"))
p2b <- ggplot(controls, aes(reorder(gsub("_", " ", control), mean_rmse), mean_rmse)) +
  geom_col(fill = palette_contract["neutral_mid"], width = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, mean_rmse - sd_rmse), ymax = mean_rmse + sd_rmse), width = 0.2, linewidth = 0.35) +
  coord_flip() + labs(x = NULL, y = "Cross-dataset RMSE (mean ± s.d.)")
variance <- fread("results/future_fate_variance_shapley.tsv")
atomic_write_tsv(variance, file.path(src_dir, "figure5_variance.tsv"))
p2c <- ggplot(variance, aes(reorder(gsub("_", " ", component), variance_fraction), variance_fraction, fill = component == "unexplained_residual")) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = unname(palette_contract["exercise"]), `TRUE` = unname(palette_contract["neutral_mid"])), guide = "none") +
  scale_y_continuous(labels = scales::percent_format()) + labs(x = NULL, y = "In-sample Shapley variance fraction")
fig2 <- p2a / p2b / p2c + plot_layout(heights=c(1.5,1,0.9)) + plot_annotation(tag_levels = "a")

# Figure 2: occupancy plasticity, entropy and model uncertainty.
dyn <- fread("results/future_unit_dynamics.tsv")
dyn[, intervention_label := label_intervention(intervention)]
atomic_write_tsv(dyn, file.path(src_dir, "figure2_dynamics.tsv"))
p3a <- ggplot(dyn, aes(reorder(intervention_label, plasticity_index_js, median), plasticity_index_js, colour = evidence_tier)) +
  geom_boxplot(outlier.shape = NA, colour = "grey50", width = 0.55, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 20261008), size = 1.1, alpha = 0.8) + coord_flip() +
  scale_colour_manual(values = c(Bronze = "#8C8C8C", Silver = "#377EB8", Gold = unname(palette_contract["gold"]))) +
  labs(x = NULL, y = "Plasticity index (JS divergence)", colour = "Evidence")
p3b <- ggplot(dyn, aes(plasticity_index_js, entropy_change_bits, colour = dataset_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2) + geom_point(size = 1.4, alpha = 0.8) +
  scale_colour_brewer(palette = "Dark2", labels=label_dataset) + labs(x = "Plasticity index (JS)", y = "Entropy change (bits)", colour = "Dataset")
unc <- fread("results/benchmark/model_disagreement_coverage_summary.tsv")
unc[, split_label := label_split(split)]
atomic_write_tsv(unc, file.path(src_dir, "figure2_uncertainty.tsv"))
p3c <- ggplot(unc, aes(reorder(split_label, mean_model_disagreement_sd), mean_model_disagreement_sd, colour = model_interval_coverage)) +
  geom_point(size = 2) + coord_flip() + scale_colour_viridis_c() +
  labs(x = NULL, y = "Model-disagreement s.d.", colour = "5–95% model\ncoverage")
fig3 <- p3a / p3b / p3c + plot_layout(heights=c(1,0.9,1.8)) + plot_annotation(tag_levels = "a")

# Figure 3: Gold genetic lineage/fusion validation.
gold <- fread("results/gold_selected_gene_results.tsv")
keep_genes <- c("Mymk", "Runx1", "Chrne", "Fbxo32", "Mtor")
gold <- gold[gene %in% keep_genes]
gold[, `:=`(
  lower = log2FoldChange - 1.96 * lfcSE,
  upper = log2FoldChange + 1.96 * lfcSE,
  contrast_label = fcase(
    grepl("^developmental", contrast), "Developmental GFP+ vs GFP−",
    grepl("newly_fused", contrast), "Newly fused GFP+ vs GFP−",
    default = "Overloaded Mymk-KO vs sham"
  )
)]
atomic_write_tsv(gold, file.path(src_dir, "figure3_gold_genes.tsv"))
p4a <- ggplot(gold, aes(log2FoldChange, gene, colour = factor(ifelse(is.na(padj), "FDR unavailable", ifelse(padj < 0.05, "FDR < 0.05", "FDR ≥ 0.05")), levels=c("FDR < 0.05", "FDR ≥ 0.05", "FDR unavailable")))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15, linewidth = 0.4, orientation = "y") + geom_point(size = 1.8) +
  facet_wrap(~ contrast_label, scales = "free_x", ncol = 1) +
  scale_colour_manual(
    values = c("FDR < 0.05" = unname(palette_contract["gold"]), "FDR ≥ 0.05" = unname(palette_contract["neutral_mid"]), "FDR unavailable" = "#222222"), drop = FALSE
  ) +
  labs(x = "log2 fold change (95% Wald CI)", y = NULL, colour = NULL)
gold_sum <- fread("results/gold_lineage_summary.tsv")
long_sum <- melt(gold_sum, id.vars = c("contrast", "biological_replicates"), measure.vars = c("n_fdr_005", "n_fdr_005_abs_lfc_1"), variable.name = "threshold", value.name = "n_genes")
long_sum[, contrast_label := fcase(
  grepl("^developmental", contrast), "Developmental lineage",
  grepl("newly_fused", contrast), "Overload newly fused",
  default = "Overloaded KO vs sham"
)]
atomic_write_tsv(long_sum, file.path(src_dir, "figure3_gold_summary.tsv"))
p4b <- ggplot(long_sum, aes(contrast_label, n_genes, fill = threshold)) +
  geom_col(position = "dodge", width = 0.7) + coord_flip() +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_fill_manual(values = c(n_fdr_005 = unname(palette_contract["neutral_mid"]), n_fdr_005_abs_lfc_1 = unname(palette_contract["gold"])), labels = c("FDR < 0.05", "FDR < 0.05, |LFC| ≥ 1")) +
  labs(x = NULL, y = "Genes (log scale)", fill = NULL)
fig4 <- p4a / p4b + plot_layout(heights = c(2, 1)) + plot_annotation(tag_levels = "a")


# Figure 4: spatial and paired-human validation.
nmj <- fread("results/nmj_section_validation.tsv")
nmj[, intervention_label := label_intervention(intervention)]
atomic_write_tsv(nmj, file.path(src_dir, "figure4_nmj_spatial.tsv"))
p5a <- ggplot(nmj, aes(intervention_label, tsc_score_local_minus_far, colour = intervention_label)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + geom_point(size = 2, position = position_jitter(width = 0.08, height = 0, seed = 20261008)) +
  scale_colour_brewer(palette = "Dark2", guide = "none") + labs(x = NULL, y = "tSC score: local minus far") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
sp <- fread("results/spatial_summary.tsv")
gene_results <- as.data.table(read.delim(xzfile("data/spatial_gene_statistics.tsv.xz"), na.strings=c("NA", "")))
gene_results[, direction_concordant_excluding_PH078 := sign(log2FoldChange) == sign(log2FoldChange_excluding_PH078)]
sig <- gene_results[!is.na(padj) & padj<0.05]
stopifnot(nrow(sig)==907L, sum(sig$direction_concordant_excluding_PH078,na.rm=TRUE)==907L,
          sum(sig$direction_concordant_excluding_PH078 & sig$padj_excluding_PH078<0.05,na.rm=TRUE)==103L)
sp_long <- data.table(stage = c("Paired FDR", "Same direction + FDR\nafter PH078 exclusion", "Exact sign-flip FDR"), genes = c(sp$n_fdr_005, sp$n_fdr_005_replicated_excluding_PH078, sp$n_exact_signflip_fdr_005))
atomic_write_tsv(sp_long, file.path(src_dir, "figure4_spatial_replication.tsv"))
p5b <- ggplot(sp_long, aes(stage, genes)) + geom_col(fill = palette_contract["exercise"], width = 0.65) +
  geom_text(aes(label = genes), vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = NULL, y = "Genes") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
human <- fread("results/human_exact_tests.tsv")
human[, cell_type_label := label_dataset(coarse_cell_type)]
atomic_write_tsv(human, file.path(src_dir, "figure4_human_paired.tsv"))
p5c <- ggplot(human, aes(difference_in_paired_change, reorder(cell_type_label, difference_in_paired_change))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) + geom_point(size = 2, colour = palette_contract["exercise"]) +
  geom_text(aes(label = sprintf("P=%.2f", exact_two_sided_p)), nudge_x = 0.006, size = 2.8, hjust = 0) +
  scale_x_continuous(limits=c(min(human$difference_in_paired_change)-0.005,max(human$difference_in_paired_change)+0.025)) +
  labs(x = "Exercise–control difference in paired change", y = NULL)
fig5 <- (p5a | p5b) / p5c + plot_annotation(tag_levels = "a")


for(name in c("Figure1_atlas_scope","Figure2_occupancy_uncertainty","Figure3_gold_lineage","Figure4_spatial_human","Figure5_generalization_controls")) {
  obj <- switch(name,Figure1_atlas_scope=fig1,Figure2_occupancy_uncertainty=fig3,Figure3_gold_lineage=fig4,Figure4_spatial_human=fig5,Figure5_generalization_controls=fig2)
  height <- switch(name,Figure1_atlas_scope=135,Figure2_occupancy_uncertainty=200,Figure3_gold_lineage=165,Figure4_spatial_human=140,Figure5_generalization_controls=165)
  ragg::agg_png(file.path("figures",paste0(name,".png")),width=183,height=height,units="mm",res=300,background="white")
  print(obj);dev.off()
}
capture.output(sessionInfo(),file="results/R_session.txt")

stopifnot(isTRUE(all.equal(ggplot_build(p5a)$data[[2]]$y,nmj$tsc_score_local_minus_far)))
stopifnot(isTRUE(all.equal(ggplot_build(p1c)$data[[2]]$y,qc$retained_fraction)))
