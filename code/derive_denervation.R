#!/usr/bin/env Rscript
# Recompute animal-level composition enumeration and summarize supplied DESeq2
# statistics. Enumeration requires exchangeable animal group labels under the
# null; it does not establish random assignment or remove confounding.
suppressPackageStartupMessages(library(data.table))
source("code/io.R")
counts <- fread("data/denervation_published_counts.tsv")
stopifnot(nrow(counts) == 78L, uniqueN(counts$sample_id) == 6L,
          uniqueN(counts$published_cell_type) == 13L,
          !anyDuplicated(counts[, .(sample_id, published_cell_type)]),
          all(is.finite(counts$n_cells)), all(counts$n_cells >= 0),
          all(counts$n_cells == round(counts$n_cells)), sum(counts$n_cells) == 29539L)
animals <- unique(counts[, .(sample_id, condition)])[order(sample_id)]
stopifnot(nrow(animals) == 6L, sum(animals$condition == "control") == 3L,
          sum(animals$condition == "denervation") == 3L,
          all(animals$sample_id[animals$condition == "control"] == paste0("GSM557158", 4:6)),
          all(animals$sample_id[animals$condition == "denervation"] == paste0("GSM557158", 7:9)))
mapping <- c(TypeI = "myonuclei", TypeIIA = "myonuclei", TypeIIx = "myonuclei", TypeIIB = "myonuclei",
             MTJ = "tendon_mtj", NMJ = "nmj_myonuclei", Satellite = "muscle_stem_myogenic",
             FAPs = "fap_mesenchymal", Pericyte = "pericyte", Adipocyte = "adipocyte",
             Endothelial = "endothelial", Macrophage = "immune", Unclassified = "other_ambiguous")
counts[, coarse_cell_type := unname(mapping[published_cell_type])]
stopifnot(!anyNA(counts$coarse_cell_type))
coarse <- counts[, .(n_cells = sum(n_cells)), by = .(sample_id, condition, coarse_cell_type)]
manifest <- fread("data/sample_manifest.tsv")[dataset_id == "denervation",
                   .(sample_id, coarse_cell_type, n_manifest = n_cells)]
check <- merge(coarse, manifest, by = c("sample_id", "coarse_cell_type"), all = TRUE)
stopifnot(nrow(check) == 60L, !anyNA(check$n_cells), !anyNA(check$n_manifest),
          all(check$n_cells == check$n_manifest), uniqueN(coarse$coarse_cell_type) == 10L)
assignments <- combn(seq_len(6L), 3L, simplify = FALSE)
assignment_table <- rbindlist(lapply(seq_along(assignments), function(j) {
  i <- assignments[[j]]
  data.table(assignment = j, denervation_samples = paste(animals$sample_id[i], collapse = ";"),
             control_samples = paste(animals$sample_id[-i], collapse = ";"))
}))
calculate <- function(x, label, representation_name) {
  x <- copy(x); setnames(x, label, "state")
  # Zero filling retains an unobserved state in every animal's denominator/test.
  grid <- CJ(sample_id = animals$sample_id, state = sort(unique(x$state)))
  z <- merge(grid, x[, .(sample_id, state, n_cells)], by = c("sample_id", "state"), all.x = TRUE)
  z[is.na(n_cells), n_cells := 0L]
  z <- merge(z, animals, by = "sample_id", all.x = TRUE)
  z[, total_cells := sum(n_cells), by = sample_id]
  z[, proportion := n_cells / total_cells]
  stopifnot(all(z$total_cells > 0), all(abs(z[, sum(proportion), by = sample_id]$V1 - 1) < 1e-12))
  tests <- list(); nulls <- list()
  for (state_name in sort(unique(z$state))) {
    a <- z[state == state_name][match(animals$sample_id, sample_id)]
    value <- a$proportion
    observed <- mean(value[a$condition == "denervation"]) - mean(value[a$condition == "control"])
    null <- vapply(assignments, function(i) mean(value[i]) - mean(value[-i]), numeric(1))
    tests[[state_name]] <- data.table(
      representation = representation_name, state = state_name,
      control_mean = mean(value[a$condition == "control"]),
      denervation_mean = mean(value[a$condition == "denervation"]),
      denervation_minus_control = observed,
      exact_two_sided_p = mean(abs(null) >= abs(observed) - 1e-15),
      n_control_animals = 3L, n_denervation_animals = 3L, n_assignments = length(null),
      minimum_attainable_two_sided_p = 2 / length(null),
      inference_condition = "animal_group_label_exchangeability_under_null")
    nulls[[state_name]] <- data.table(representation = representation_name, state = state_name,
                                     assignment = seq_along(null), mean_difference = null)
  }
  result <- rbindlist(tests)
  result[, fdr_bh := p.adjust(exact_two_sided_p, method = "BH")]
  z[, representation := representation_name]
  list(composition = z, tests = result, nulls = rbindlist(nulls))
}
published <- calculate(counts, "published_cell_type", "published_13_labels")
coarsened <- calculate(coarse, "coarse_cell_type", "observed_coarse_10_categories")
# schwann_glial and smooth_muscle are in the global 12-category dictionary but
# cannot be recovered separately from these published identities. They are not
# added as two identically zero hypothesis tests to this cohort's analysis.
atomic_write_tsv(rbindlist(list(published$composition, coarsened$composition)), "results/denervation_animal_composition.tsv")
atomic_write_tsv(rbindlist(list(published$tests, coarsened$tests)), "results/denervation_composition_exact_tests.tsv")
atomic_write_tsv(assignment_table, "results/denervation_composition_assignments.tsv")
atomic_write_tsv(rbindlist(list(published$nulls, coarsened$nulls)), "results/denervation_composition_null_statistics.tsv")
atomic_write_tsv(data.table(published_cell_type = names(mapping), coarse_cell_type = unname(mapping)), "results/denervation_label_mapping.tsv")
# Animal-pseudobulk estimates are upstream fitted inputs. Here only count
# summaries are recomputed; adjusted P values and coefficients are not refitted.
de <- as.data.table(read.delim(xzfile("data/denervation_gene_statistics.tsv.xz"), na.strings = c("NA", ""), check.names = FALSE))
stopifnot(!anyDuplicated(de[, .(cell_type, gene)]), uniqueN(de$cell_type) == 7L,
          all(is.na(de$padj) | (de$padj >= 0 & de$padj <= 1)))
summary <- de[, .(
  n_gene_rows = .N, n_nonmissing_padj = sum(!is.na(padj)),
  n_fdr_005 = sum(padj < 0.05, na.rm = TRUE),
  n_fdr_005_abs_lfc_1 = sum(padj < 0.05 & abs(log2FoldChange) >= 1, na.rm = TRUE)
), by = cell_type][order(cell_type)]
summary[, `:=`(contrast = "denervation_vs_control", n_control_animals = 3L,
               n_denervation_animals = 3L, inference_unit = "animal_pseudobulk",
               statistics_origin = "supplied_previously_fitted_DESeq2_statistics",
               reported_named_type_in_manuscript = cell_type != "Unclassified")]
atomic_write_tsv(summary, "results/denervation_pseudobulk_summary.tsv")
print(published$tests)
print(coarsened$tests)
print(summary)
