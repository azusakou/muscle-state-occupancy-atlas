#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(glmnet)
  library(nnet)
  library(xgboost)
})
source("code/io.R")
source("code/train_only_encoding.R")

set.seed(20260715)
tasks <- as.data.table(read.delim(gzfile("data/future_tasks.tsv.gz"), na.strings = c("NA", ""), check.names=FALSE))[primary_task == TRUE]
if (!nrow(tasks)) stop("No primary future tasks")
tasks[, unit_id := paste(dataset_id, biological_unit_id, sep = "||")]
tasks[, composition_id := paste(dataset_id, sample_id, sep = "||")]
tasks[, target := pmin(pmax(proportion, 1e-5), 1 - 1e-5)]
tasks[, target_logit := qlogis(target)]
tasks[, baseline_logit := qlogis(pmin(pmax(baseline_proportion, 1e-5), 1 - 1e-5))]
if (!all(c("time_days", "age_months") %in% names(tasks))) stop("Canonical time_days/age_months fields are missing; rerun scripts/42_build_future_tasks.R")
tasks[, log_time := log1p(time_days)]
tasks[, time_bin := factor(ifelse(time_days <= 3, "early", "late"))]
tasks[, age_numeric := age_months]
tasks[, age_bin := factor(ifelse(is.na(age_months), "unknown", ifelse(age_months <= 12, "younger", "older")))]
for (field in c("coarse_cell_type", "intervention", "species", "dataset_id")) tasks[, (field) := factor(get(field))]


splits <- list()
add_split <- function(name, test) {
  test <- which(test)
  train <- setdiff(seq_len(nrow(tasks)), test)
  if (length(test) && length(train) && uniqueN(tasks$unit_id[test]) >= 1L && uniqueN(tasks$unit_id[train]) >= 2L)
    splits[[name]] <<- list(train = train, test = test)
}
for (value in levels(tasks$dataset_id)) add_split(paste0("leave_dataset__", value), tasks$dataset_id == value)
for (value in levels(tasks$intervention)) add_split(paste0("leave_intervention__", value), tasks$intervention == value)
for (value in levels(tasks$species)) add_split(paste0("leave_species__", gsub(" ", "_", value)), tasks$species == value)
add_split("time_extrapolation__gt3_days", tasks$time_days > 3)
for (value in levels(tasks$age_bin)) add_split(paste0("leave_age_bin__", value), tasks$age_bin == value)

fallback_group_mean <- function(train, test, fields) {
  global <- mean(tasks$target_logit[train])
  lookup <- tasks[train, .(prediction = mean(target_logit)), by = fields]
  query <- copy(tasks[test])
  query[, row_order__ := .I]
  query <- merge(query, lookup, by = fields, all.x = TRUE, sort = FALSE)
  setorder(query, row_order__)
  prediction <- query$prediction
  if (anyNA(prediction) && length(fields) > 1L) {
    cell_lookup <- tasks[train, .(cell_fallback = mean(target_logit)), by = coarse_cell_type]
    prediction[is.na(prediction)] <- cell_lookup$cell_fallback[
      match(tasks$coarse_cell_type[test][is.na(prediction)], cell_lookup$coarse_cell_type)
    ]
  }
  prediction[is.na(prediction)] <- global
  prediction
}

ols_predict <- function(X, y, train, test, ridge = 1e-7) {
  xtx <- crossprod(X[train, , drop = FALSE])
  penalty <- diag(ridge, ncol(X)); penalty[1, 1] <- 0
  beta <- solve(xtx + penalty, crossprod(X[train, , drop = FALSE], y[train]))
  drop(X[test, , drop = FALSE] %*% beta)
}

glmnet_predict <- function(X, y, train, test, alpha, lambda) {
  fit <- glmnet(X[train, , drop = FALSE], y[train], alpha = alpha, lambda = lambda,
                standardize = TRUE, intercept = TRUE)
  drop(predict(fit, X[test, , drop = FALSE], s = lambda))
}

knn_predict <- function(X, y, train, test, k = 5L) {
  center <- colMeans(X[train, , drop = FALSE])
  scale <- apply(X[train, , drop = FALSE], 2, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  tr <- sweep(sweep(X[train, , drop = FALSE], 2, center), 2, scale, "/")
  te <- sweep(sweep(X[test, , drop = FALSE], 2, center), 2, scale, "/")
  vapply(seq_len(nrow(te)), function(i) {
    distance <- rowSums((tr - matrix(te[i, ], nrow(tr), ncol(tr), byrow = TRUE))^2)
    mean(y[train][order(distance)[seq_len(min(k, length(distance)))]] )
  }, numeric(1))
}

model_names <- c(
  "persistence", "global_mean", "celltype_mean", "intervention_celltype_mean",
  "timebin_celltype_mean", "species_celltype_mean", "agebin_celltype_mean",
  "ols_baseline", "ols_context", "ols_additive", "ols_interactions",
  "ridge_lambda_0.01", "ridge_lambda_0.1", "ridge_lambda_1", "ridge_lambda_10",
  "lasso_lambda_0.1", "elasticnet_lambda_0.1", "knn_k5", "neuralnet_hidden5", "xgboost_shallow"
)

prediction_tables <- list()
preprocessing_audit <- list()
convergence_audit <- list()
counter <- 0L
for (split_name in names(splits)) {
  split <- splits[[split_name]]
  train <- split$train; test <- split$test
  stopifnot(!length(intersect(tasks$unit_id[train],tasks$unit_id[test])))
  f <- fold_prepare(copy(tasks[train]),copy(tasks[test]))
  tasks$age_numeric[train] <- f$train$age_numeric; tasks$age_numeric[test] <- f$test$age_numeric
  build_x <- function(formula) {
    z <- fold_matrix(f$train,f$test,formula)
    z[order(c(train,test)),,drop=FALSE]
  }
  X_base <- build_x(~ baseline_logit)
  X_context <- build_x(~ coarse_cell_type + intervention + log_time + age_numeric + species)
  X_add <- build_x(~ baseline_logit + coarse_cell_type + intervention + log_time + age_numeric + species)
  X_interaction <- build_x(~ baseline_logit * coarse_cell_type + intervention * coarse_cell_type + log_time * coarse_cell_type + age_numeric + species)
  preprocessing_audit[[split_name]] <- data.table(split=split_name,training_age_median=f$median,n_training_units=uniqueN(tasks$unit_id[train]),n_test_units=uniqueN(tasks$unit_id[test]))
  y <- tasks$target_logit
  predictions <- list(
    persistence = tasks$baseline_logit[test],
    global_mean = rep(mean(y[train]), length(test)),
    celltype_mean = fallback_group_mean(train, test, "coarse_cell_type"),
    intervention_celltype_mean = fallback_group_mean(train, test, c("intervention", "coarse_cell_type")),
    timebin_celltype_mean = fallback_group_mean(train, test, c("time_bin", "coarse_cell_type")),
    species_celltype_mean = fallback_group_mean(train, test, c("species", "coarse_cell_type")),
    agebin_celltype_mean = fallback_group_mean(train, test, c("age_bin", "coarse_cell_type")),
    ols_baseline = ols_predict(X_base, y, train, test),
    ols_context = ols_predict(X_context, y, train, test),
    ols_additive = ols_predict(X_add, y, train, test),
    ols_interactions = ols_predict(X_interaction, y, train, test),
    ridge_lambda_0.01 = glmnet_predict(X_interaction, y, train, test, 0, 0.01),
    ridge_lambda_0.1 = glmnet_predict(X_interaction, y, train, test, 0, 0.1),
    ridge_lambda_1 = glmnet_predict(X_interaction, y, train, test, 0, 1),
    ridge_lambda_10 = glmnet_predict(X_interaction, y, train, test, 0, 10),
    lasso_lambda_0.1 = glmnet_predict(X_interaction, y, train, test, 1, 0.1),
    elasticnet_lambda_0.1 = glmnet_predict(X_interaction, y, train, test, 0.5, 0.1),
    knn_k5 = knn_predict(X_add, y, train, test, 5)
  )
  set.seed(20260715)
  neural <- nnet(
    x = X_add[train, -1, drop = FALSE], y = y[train], size = 5,
    linout = TRUE, decay = 0.01, maxit = 5000, trace = FALSE
  )
  convergence_audit[[split_name]] <- data.table(split=split_name,model="neuralnet_hidden5",convergence=neural$convergence,objective=neural$value)
  predictions$neuralnet_hidden5 <- drop(predict(neural, X_add[test, -1, drop = FALSE], type = "raw"))
  booster <- xgb.train(
    params = list(
      objective = "reg:squarederror", max_depth = 2, eta = 0.05,
      subsample = 0.8, colsample_bytree = 0.8, nthread = 1, seed = 20260715
    ),
    data = xgb.DMatrix(X_add[train, -1, drop = FALSE], label = y[train]),
    nrounds = 100, verbose = 0
  )
  predictions$xgboost_shallow <- predict(booster, X_add[test, -1, drop = FALSE])

  stopifnot(all(vapply(predictions,function(x) length(x)==length(test)&&all(is.finite(x)),logical(1))))
  for (model in model_names) {
    counter <- counter + 1L
    table <- copy(tasks[test])
    table[, `:=`(
      split = split_name,
      model = model,
      prediction_raw = plogis(pmin(pmax(predictions[[model]], -20), 20))
    )]
    table[, prediction := prediction_raw / sum(prediction_raw), by = composition_id]
    prediction_tables[[counter]] <- table
  }
}

pred <- rbindlist(prediction_tables, fill = TRUE)
pred[, squared_error := (prediction - proportion)^2]
pred[, absolute_error := abs(prediction - proportion)]
metrics <- pred[, .(
  n_test_units = uniqueN(unit_id),
  n_test_rows = .N,
  mae = mean(absolute_error),
  rmse = sqrt(mean(squared_error)),
  spearman = suppressWarnings(cor(prediction, proportion, method = "spearman"))
), by = .(split, model)]

js <- pred[, {
  p <- pmax(proportion, 1e-12); p <- p / sum(p)
  q <- pmax(prediction, 1e-12); q <- q / sum(q)
  m <- (p + q) / 2
  .(
    js_divergence = 0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m)),
    dominant_type_correct = coarse_cell_type[which.max(p)] == coarse_cell_type[which.max(q)]
  )
}, by = .(split, model, composition_id)]
js_summary <- js[, .(
  mean_js_divergence = mean(js_divergence),
  median_js_divergence = median(js_divergence),
  dominant_type_accuracy = mean(dominant_type_correct)
), by = .(split, model)]
metrics <- merge(metrics, js_summary, by = c("split", "model"))
setorder(metrics, split, mean_js_divergence, mae)

out <- "results/benchmark"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
atomic_write_tsv(rbindlist(convergence_audit),file.path(out,"model_convergence_audit.tsv"))
atomic_write_tsv(rbindlist(preprocessing_audit),file.path(out,"preprocessing_audit.tsv"))
atomic_write_tsv(pred, file.path(out, "benchmark_predictions.tsv.gz"))
atomic_write_tsv(metrics, file.path(out, "benchmark_metrics.tsv"))
atomic_write_tsv(js, file.path(out, "benchmark_unit_js.tsv"))
print(metrics[, head(.SD, 3), by = split])
