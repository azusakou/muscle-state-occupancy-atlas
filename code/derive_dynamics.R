#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source("code/io.R")

tasks <- as.data.table(read.delim(gzfile("data/future_tasks.tsv.gz"),na.strings=c("NA","")))[primary_task %in% c(TRUE,"True")]
if (!nrow(tasks)) stop("No primary future tasks")

# These are state-occupancy changes between independent or linked samples. They
# are deliberately not called single-cell transition probabilities unless the
# registry provides lineage/longitudinal evidence.
tasks[, baseline_p := pmax(baseline_proportion, 1e-12)]
tasks[, future_p := pmax(proportion, 1e-12)]
tasks[, occupancy_change := proportion - baseline_proportion]
tasks[, log2_occupancy_ratio := log2(future_p / baseline_p)]
tasks[, state_expansion_probability_proxy := pmax(occupancy_change, 0)]
tasks[, state_extinction_probability_proxy := pmax(-occupancy_change, 0)]

unit_summary <- tasks[, {
  p <- baseline_p / sum(baseline_p)
  q <- future_p / sum(future_p)
  m <- (p + q) / 2
  entropy <- function(x) -sum(x * log2(x))
  js <- 0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
  .(
    baseline_entropy_bits = entropy(p),
    future_entropy_bits = entropy(q),
    entropy_change_bits = entropy(q) - entropy(p),
    plasticity_index_js = js,
    total_variation_distance = 0.5 * sum(abs(q - p)),
    dominant_baseline_state = coarse_cell_type[which.max(p)],
    dominant_future_state = coarse_cell_type[which.max(q)],
    dominant_state_changed = coarse_cell_type[which.max(p)] != coarse_cell_type[which.max(q)]
  )
}, by = .(
  dataset_id, biological_unit_id, sample_id, species, intervention,
  time_value, time_unit, evidence_tier, temporal_link_type, claim_scope
)]


atomic_write_tsv(unit_summary,"results/future_unit_dynamics.tsv")
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
source("code/io.R")
set.seed(20260715)

x <- as.data.table(read.delim(gzfile("data/future_tasks.tsv.gz"),na.strings=c("NA","")))[primary_task %in% c(TRUE,"True")]
x[, target_logit := qlogis(pmin(pmax(proportion, 1e-5), 1 - 1e-5))]
x[, baseline_logit := qlogis(pmin(pmax(baseline_proportion, 1e-5), 1 - 1e-5))]
if (!all(c("time_days", "age_months") %in% names(x))) stop("Canonical time_days/age_months fields are missing; rerun scripts/42_build_future_tasks.R")
x[, log_time := log1p(time_days)]
x[, age_numeric := fifelse(is.na(age_months), median(age_months, na.rm = TRUE), age_months)]

groups <- list(
  baseline_composition = "baseline_logit",
  cell_type_identity = "coarse_cell_type",
  intervention_specific = "intervention",
  time_dependent = "log_time",
  species_age_context = "species + age_numeric",
  study_context = "dataset_id"
)
r2_for <- function(selected) {
  if (!length(selected)) return(0)
  terms <- unlist(groups[selected], use.names = FALSE)
  fit <- lm(as.formula(paste("target_logit ~", paste(terms, collapse = " + "))), data = x)
  summary(fit)$r.squared
}

# Exact Shapley decomposition over all 6! group orders removes arbitrary
# dependence on the order in which covariate blocks enter the model.
permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
  do.call(rbind, lapply(seq_along(v), function(i) cbind(v[i], permutations(v[-i]))))
}
orders <- permutations(names(groups))
contrib <- setNames(numeric(length(groups)), names(groups))
for (i in seq_len(nrow(orders))) {
  selected <- character()
  previous <- 0
  for (g in orders[i, ]) {
    current <- r2_for(c(selected, g))
    contrib[g] <- contrib[g] + current - previous
    selected <- c(selected, g)
    previous <- current
  }
}
contrib <- contrib / nrow(orders)
full_r2 <- r2_for(names(groups))
decomposition <- data.table(
  component = c(names(contrib), "unexplained_residual"),
  variance_fraction = c(unname(contrib), 1 - full_r2),
  estimand = c(
    "baseline compositional state (not cell-intrinsic expression)", "harmonized cell-type identity",
    "intervention label", "sampling time", "species and age", "dataset/study context",
    "unexplained variation; not identifiable as biological stochasticity"
  )
)
decomposition[, fraction_of_explained := fifelse(component == "unexplained_residual", NA_real_, variance_fraction / full_r2)]


atomic_write_tsv(decomposition,"results/future_fate_variance_shapley.tsv")
