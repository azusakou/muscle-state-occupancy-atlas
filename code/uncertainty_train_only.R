#!/usr/bin/env Rscript


suppressPackageStartupMessages(library(data.table))
source("code/io.R")
set.seed(20260715)

pred <- as.data.table(read.delim(gzfile("results/benchmark/benchmark_predictions.tsv.gz"),na.strings=c("NA","")))
if (!nrow(pred)) stop("No benchmark predictions")

# Descriptive model-specification sensitivity: disagreement across the deliberately diverse specification set for
# the same held-out analysis unit and cell state. This is not relabelled as
# aleatoric uncertainty.
ensemble <- pred[, .(
  observed = unique(proportion)[1],
  prediction_mean = mean(prediction),
  prediction_sd_across_models = sd(prediction),
  prediction_q05 = quantile(prediction, 0.05),
  prediction_q50 = quantile(prediction, 0.50),
  prediction_q95 = quantile(prediction, 0.95),
  n_models = uniqueN(model)
), by = .(split, unit_id, coarse_cell_type)]
ensemble[, covered_by_model_q05_q95 := observed >= prediction_q05 & observed <= prediction_q95]

# Conditional test-unit bootstrap, holding fitted models fixed, preserves all compositional rows from a held-out unit.
metric_boot <- list(); k <- 0L
for (s in unique(pred$split)) for (m in unique(pred[split == s]$model)) {
  x <- pred[split == s & model == m]
  units <- unique(x$unit_id)
  if (length(units) < 2L) next
  unit_loss <- x[, .(
    n_rows = .N,
    absolute_error_sum = sum(abs(prediction - proportion)),
    squared_error_sum = sum((prediction - proportion)^2)
  ), by = unit_id]
  for (b in seq_len(250L)) {
    sampled <- unit_loss[sample.int(.N, .N, replace = TRUE)]
    k <- k + 1L
    metric_boot[[k]] <- data.table(
      split = s, model = m, bootstrap = b,
      mae = sum(sampled$absolute_error_sum) / sum(sampled$n_rows),
      rmse = sqrt(sum(sampled$squared_error_sum) / sum(sampled$n_rows))
    )
  }
}
boot <- rbindlist(metric_boot)
intervals <- boot[, .(
  mae_q025 = quantile(mae, 0.025), mae_median = median(mae), mae_q975 = quantile(mae, 0.975),
  rmse_q025 = quantile(rmse, 0.025), rmse_median = median(rmse), rmse_q975 = quantile(rmse, 0.975)
), by = .(split, model)]

coverage <- ensemble[, .(
  n_units = uniqueN(unit_id), n_rows = .N,
  mean_model_disagreement_sd = mean(prediction_sd_across_models),
  model_interval_coverage = mean(covered_by_model_q05_q95)
), by = split]

out <- "results/benchmark"
atomic_write_tsv(ensemble, file.path(out, "model_disagreement_uncertainty.tsv.gz"))
atomic_write_tsv(intervals, file.path(out, "biological_unit_bootstrap_metric_intervals.tsv"))
atomic_write_tsv(coverage, file.path(out, "model_disagreement_coverage_summary.tsv"))
print(coverage)
