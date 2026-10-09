suppressPackageStartupMessages(library(data.table))
source("code/io.R")
set.seed(20260715)
# Lineage estimates are refitted from counts; spatial estimates are processed inputs.
gold <- fread("results/lineage_gene_statistics.tsv")
genes <- c("Mymk","Runx1","Chrne","Fbxo32","Mtor")
stopifnot(all(genes %in% gold$gene))
atomic_write_tsv(gold[gene %in% genes],"results/gold_selected_gene_results.tsv")
summary <- gold[, .(n_genes_tested=.N,n_fdr_005=sum(padj<.05,na.rm=TRUE),n_fdr_005_abs_lfc_1=sum(padj<.05 & abs(log2FoldChange)>=1,na.rm=TRUE)),by=contrast]
summary$biological_replicates <- ifelse(grepl("developmental|newly_fused",summary$contrast),"3 paired pooled preparations; preparation-adjusted","3 overloaded KO vs 4 non-overloaded sham preparations; batch-adjusted; only one sham in KO batch")
atomic_write_tsv(summary,"results/gold_lineage_summary.tsv")
sp <- as.data.table(read.delim(xzfile("data/spatial_gene_statistics.tsv.xz"), na.strings=c("NA", "")))
same <- sign(sp$log2FoldChange)==sign(sp$log2FoldChange_excluding_PH078)
atomic_write_tsv(data.table(n_fdr_005=sum(sp$padj<.05,na.rm=TRUE),
 n_direction_retained=sum(sp$padj<.05 & same,na.rm=TRUE),
 n_fdr_005_replicated_excluding_PH078=sum(sp$padj<.05 & sp$padj_excluding_PH078<.05 & same,na.rm=TRUE),
 n_exact_signflip_fdr_005=sum(sp$exact_signflip_padj<.05,na.rm=TRUE)),"results/spatial_summary.tsv")
# Recompute all 15 participant assignments from observed paired changes.
wide <- fread("data/human_paired_changes.tsv")
subjects <- unique(wide[,.(subject,Group)])
assignments <- combn(subjects$subject,sum(subjects$Group=="Exercise"),simplify=FALSE)
tests <- wide[, {
 observed <- mean(post_minus_pre[Group=="Exercise"])-mean(post_minus_pre[Group=="Control"])
 null <- vapply(assignments,function(ids) mean(post_minus_pre[subject %in% ids])-mean(post_minus_pre[!subject %in% ids]),numeric(1))
 .(exercise_mean_paired_change=mean(post_minus_pre[Group=="Exercise"]),
 control_mean_paired_change=mean(post_minus_pre[Group=="Control"]),
 difference_in_paired_change=observed,exact_two_sided_p=mean(abs(null)>=abs(observed)-1e-15),
 n_exercise_subjects=sum(Group=="Exercise"),n_control_subjects=sum(Group=="Control"),n_exact_assignments=length(null))
},by=coarse_cell_type]
tests$fdr_bh <- p.adjust(tests$exact_two_sided_p,"BH")
atomic_write_tsv(tests,"results/human_exact_tests.tsv")
# Recompute the random-anchor diagnostic from processed spot scores and coordinates.
spots <- fread("data/nmj_spots.tsv")
parts <- list()
for(sample in unique(spots$sample_id)) {
 z <- spots[sample_id==sample]
 xy <- as.matrix(z[,.(array_row,array_col)])
 distance <- as.matrix(dist(xy))
 neighbors <- lapply(seq_len(nrow(z)),function(i) which(distance[,i]<=4))
 statistic <- function(anchors) {
   local <- rep(FALSE,nrow(z))
   for(i in anchors)local[neighbors[[i]]]<-TRUE
   mean(z$terminal_schwann_score[local])-mean(z$terminal_schwann_score[!local])
 }
 observed <- statistic(which(z$chrne_positive))
 null <- replicate(500,statistic(sample.int(nrow(z),sum(z$chrne_positive),replace=FALSE)))
 parts[[sample]] <- data.table(sample_id=sample,intervention=z$intervention[1],
 tsc_score_local_minus_far=observed,spatial_label_permutation_p=(1+sum(null>=observed))/501,
 permutation_n=500)
}
atomic_write_tsv(rbindlist(parts),"results/nmj_section_validation.tsv")
