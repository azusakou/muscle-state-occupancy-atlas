#!/usr/bin/env Rscript
# Refit three sample-level contrasts from gene-collapsed integer counts.
# GFP+/- fractions share three pooled preparations within each lineage experiment.
# KO vs SHAM contrasts loading/genotype contexts, not an isolated fusion effect.
suppressPackageStartupMessages({library(DESeq2); library(data.table)})
source("code/io.R")
x <- readRDS("data/lineage_counts.rds")
m <- x$sample_metadata
stopifnot(identical(colnames(x$collapsed_counts), as.character(m$sample)),
          !anyDuplicated(m$sample), !anyDuplicated(rownames(x$collapsed_counts)),
          all(x$collapsed_counts >= 0))
definitions <- list(
  list("DEV_POS", "DEV_NEG", "developmental_GFP_positive_vs_negative", TRUE),
  list("MOV_POS", "MOV_NEG", "overload_newly_fused_GFP_positive_vs_preexisting_GFP_negative", TRUE),
  list("MYMK_KO", "SHAM", "overload_fusion_blockade_vs_sham_batch_adjusted", FALSE))
outputs <- audits <- list()
for (d in definitions) {
  meta <- droplevels(m[m$condition %in% d[1:2], , drop=FALSE])
  meta$condition <- relevel(factor(meta$condition), d[[2]])
  meta$preparation_pair <- factor(meta$preparation_pair)
  meta$batch <- factor(meta$batch)
  design <- if(d[[4]]) ~preparation_pair + condition else ~batch + condition
  if(d[[4]]) stopifnot(all(table(meta$preparation_pair, meta$condition)==1L))
  mat <- model.matrix(design, meta)
  stopifnot(qr(mat)$rank==ncol(mat))
  counts <- x$collapsed_counts[, meta$sample, drop=FALSE]
  keep <- rowSums(counts >= 10) >= 2
  dds <- DESeqDataSetFromMatrix(counts[keep,,drop=FALSE], meta, design=design)
  dds <- DESeq(dds, quiet=TRUE)
  # If failures occur, rerun coefficient fitting with a larger optimization budget.
  initial_failed <- sum(!mcols(dds)$betaConv, na.rm=TRUE)
  if(initial_failed) dds <- nbinomWaldTest(dds, maxit=1000, quiet=TRUE)
  stopifnot(all(mcols(dds)$betaConv))
  z <- as.data.table(as.data.frame(results(dds, contrast=c("condition",d[[1]],d[[2]]), alpha=0.05)),keep.rownames="gene")
  z[, contrast := d[[3]]]
  outputs[[d[[3]]]] <- z
  audits[[d[[3]]]] <- data.table(contrast=d[[3]],design=paste(deparse(design),collapse=""),
    n_samples=nrow(meta),n_parameters=ncol(mat),design_rank=qr(mat)$rank,
    n_genes=nrow(z),initial_failed_coefficients=initial_failed,
    final_failed_coefficients=sum(!mcols(dds)$betaConv),
    n_outlier_replacements=sum(mcols(dds)$replace,na.rm=TRUE),
    independent_filtering_alpha=0.05,
    deseq2_version=as.character(packageVersion("DESeq2")))
}
atomic_write_tsv(rbindlist(outputs),"results/lineage_gene_statistics.tsv")
atomic_write_tsv(rbindlist(audits),"results/lineage_design_audit.tsv")
capture.output(sessionInfo(),file="results/R_session_lineage.txt")
