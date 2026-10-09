#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})
source("code/io.R")
source("code/train_only_encoding.R")
set.seed(20260715)

exclude_ambiguous <- "--exclude-ambiguous" %in% commandArgs(trailingOnly=TRUE)
input <- if(exclude_ambiguous) "results/tasks_exclude_ambiguous.tsv.gz" else "data/future_tasks.tsv.gz"
tasks <- as.data.table(read.delim(gzfile(input), na.strings = c("NA", ""), check.names=FALSE))[primary_task %in% c(TRUE,"True")]
if (!all(c("time_days", "age_months") %in% names(tasks))) stop("Canonical time_days/age_months fields are missing; rerun scripts/42_build_future_tasks.R")
tasks[, unit_id := paste(dataset_id, biological_unit_id, sep = "||")]
tasks[, composition_id := paste(dataset_id, sample_id, sep = "||")]
tasks[, target := qlogis(pmin(pmax(proportion, 1e-5), 1 - 1e-5))]
tasks[, baseline_logit := qlogis(pmin(pmax(baseline_proportion, 1e-5), 1 - 1e-5))]
tasks[, log_time := log1p(time_days)]
tasks[, age_numeric := age_months]
tasks[, log_depth := log1p(if ("capture_total_cells" %in% names(tasks)) capture_total_cells else total_cells)]
for (v in c("coarse_cell_type", "intervention", "species", "dataset_id")) tasks[, (v) := factor(get(v))]

make_x <- function(d, formula) model.matrix(formula, d)
predict_ridge <- function(train, test, formula) {
  f <- fold_prepare(train,test)
  x <- fold_matrix(f$train,f$test,formula)
  tr <- seq_len(nrow(train)); te <- nrow(train) + seq_len(nrow(test))
  fit <- glmnet(x[tr, , drop = FALSE], train$target, alpha = 0, lambda = 0.1, standardize = TRUE)
  drop(predict(fit, x[te, , drop = FALSE], s = 0.1))
}
score <- function(test, pred) {
  scored <- copy(test)
  scored[, prediction_raw := plogis(pmin(pmax(pred, -20), 20))]
  scored[, prediction := prediction_raw / sum(prediction_raw), by = composition_id]
  sqrt(mean((scored$prediction - scored$proportion)^2))
}

full_formula <- ~ baseline_logit + coarse_cell_type + intervention + log_time + age_numeric + species

results <- list(); k <- 0L
for (held in levels(tasks$dataset_id)) {
  train <- copy(tasks[dataset_id != held]); test <- copy(tasks[dataset_id == held])
  if (!nrow(train) || !nrow(test)) next
  stopifnot(!length(intersect(train$unit_id,test$unit_id)))
  controls <- list(
    observed_baseline_context = full_formula,
    cell_composition_only = ~ baseline_logit + coarse_cell_type,
    future_capture_cell_count_only = ~ log_depth + coarse_cell_type,
    study_identity_only = ~ dataset_id + coarse_cell_type,
    species_identity_only = ~ species + coarse_cell_type
  )
  for (nm in names(controls)) {
    k <- k + 1L
    results[[k]] <- data.table(
      held_out_dataset = held, control = nm, replicate = 0L,
      rmse = score(test, predict_ridge(train, test, controls[[nm]]))
    )
  }
  for (replicate in 1:25) {
    perm <- copy(train)
    # Permute complete baseline-composition vectors between analysis units so
    # each donor vector still sums to one and retains its cell-state covariance.
    unit_map <- unique(perm[, .(unit_id)])
    unit_map[, donor_unit_id := sample(unit_id)]
    donor <- unique(perm[, .(donor_unit_id = unit_id, coarse_cell_type, donor_baseline_logit = baseline_logit)])
    perm <- merge(perm, unit_map, by = "unit_id", all.x = TRUE, sort = FALSE)
    perm <- merge(perm, donor, by = c("donor_unit_id", "coarse_cell_type"), all.x = TRUE, sort = FALSE)
    if (anyNA(perm$donor_baseline_logit)) stop("Whole-vector baseline permutation produced missing states")
    perm[, baseline_logit := donor_baseline_logit]
    k <- k + 1L
    results[[k]] <- data.table(
      held_out_dataset = held, control = "baseline_vector_permutation", replicate = replicate,
      rmse = score(test, predict_ridge(perm, test, full_formula))
    )
    perm <- copy(train)
    unit_map <- unique(perm[, .(unit_id, intervention)])
    unit_map[, intervention := sample(intervention)]
    perm[, intervention := NULL]
    perm <- merge(perm, unit_map, by = "unit_id", all.x = TRUE, sort = FALSE)
    perm[, intervention := factor(intervention, levels = levels(tasks$intervention))]
    k <- k + 1L
    results[[k]] <- data.table(
      held_out_dataset = held, control = "intervention_label_permutation", replicate = replicate,
      rmse = score(test, predict_ridge(perm, test, full_formula))
    )
    perm <- copy(train)
    unit_map <- unique(perm[, .(unit_id, log_time)])
    unit_map[, log_time := sample(log_time)]
    perm[, log_time := NULL]
    perm <- merge(perm, unit_map, by = "unit_id", all.x = TRUE, sort = FALSE)
    k <- k + 1L
    results[[k]] <- data.table(
      held_out_dataset = held, control = "timepoint_permutation", replicate = replicate,
      rmse = score(test, predict_ridge(perm, test, full_formula))
    )
  }
}
numeric_results <- rbindlist(results, fill = TRUE)
# Dataset-paired comparisons use one mean RMSE per held-out dataset and control.
# Overlapping training sets preclude treating fold errors as independent.
split_means <- numeric_results[, .(rmse = mean(rmse)), by = .(held_out_dataset, control)]
replicate_counts <- numeric_results[, .(n_replicates = .N), by = control]
summary <- split_means[, .(
  n_splits = .N, mean_rmse = mean(rmse), sd_rmse = sd(rmse), median_rmse = median(rmse)
), by = control]
summary <- merge(summary, replicate_counts, by = "control", sort = FALSE)[order(mean_rmse)]
observed <- split_means[control == "observed_baseline_context", .(held_out_dataset, observed_rmse = rmse)]
paired <- merge(split_means[control != "observed_baseline_context"], observed, by = "held_out_dataset")
# Training sets overlap: paired differences are descriptive, without fold-based P values.
paired[, rmse_advantage_of_observed := rmse - observed_rmse]
paired_summary <- paired[, .(
  n_datasets = .N,
  mean_rmse_advantage_of_observed = mean(rmse_advantage_of_observed),
  observed_better_count = sum(observed_rmse < rmse),
  minimum_advantage = min(rmse_advantage_of_observed),
  maximum_advantage = max(rmse_advantage_of_observed)
), by = control][order(-mean_rmse_advantage_of_observed)]

audit <- data.table(
  requested_control = c(
    "baseline composition-vector permutation", "intervention label permutation", "timepoint permutation",
    "clone label permutation", "random gene set", "pathway-size-matched gene set",
    "cell-count-only model", "sequencing-depth-only model", "study-identity-only model",
    "batch-only model", "post-treatment leakage control", "pseudo-time leakage control",
    "patient-identity leakage control", "muscle-group confounding control",
    "species-identity confounding control", "exercise-dose shuffling",
    "random optimal-transport coupling", "random spatial-neighborhood graph"
  ),
  status = c(
    "performed", "performed", "performed", "not_applicable_no_clone_labels",
    "not_applicable_no_gene_level_predictor", "not_applicable_no_gene_level_predictor",
    "performed_future_capture_cell_count_nuisance_diagnostic", "not_performed_no_read_or_UMI_depth_predictor", "performed",
    "not_estimable_batch_not_harmonized", "observed_model_baseline_and_design_only_depth_control_uses_future_capture",
    "passed_no_pseudotime_predictor", "passed_split_by_declared_analysis_unit",
    "not_estimable_muscle_constant_or_missing_within_study", "performed",
    "not_estimable_no_continuous_dose", "not_applicable_no_OT_claim",
    "not_applicable_spatial_test_uses_observed_coordinates_and_subject_sign_flips"
  )
)

out <- if(exclude_ambiguous) "results/controls_exclude_ambiguous" else "results/benchmark"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
atomic_write_tsv(numeric_results, file.path(out, "negative_control_split_results.tsv"))
atomic_write_tsv(summary, file.path(out, "negative_control_summary.tsv"))
atomic_write_tsv(paired, file.path(out, "negative_control_paired_dataset_differences.tsv"))
atomic_write_tsv(paired_summary, file.path(out, "negative_control_paired_comparisons.tsv"))
atomic_write_tsv(audit, file.path(out, "leakage_confounding_audit.tsv"))
print(summary)
print(paired_summary)
print(audit)
