#!/usr/bin/env Rscript
# R1.5: design identifiability and descriptive within-batch influence diagnostics.
# No new P values are assigned to the three-versus-one restricted comparison.
x <- readRDS('data/lineage_counts.rds')
out <- 'results'
m <- x$sample_metadata
stopifnot(identical(colnames(x$collapsed_counts),as.character(m$sample)),!anyDuplicated(m$sample),!anyDuplicated(rownames(x$collapsed_counts)))
keep <- m$condition %in% c('MYMK_KO','SHAM'); m <- droplevels(m[keep,]); counts <- round(x$collapsed_counts[,keep])
m$condition <- relevel(factor(m$condition),'SHAM');m$batch <- factor(m$batch)
write.table(as.data.frame(table(m$batch,m$condition)),file.path(out,'fusion_batch_design.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
X <- model.matrix(~batch+condition,m)
without_anchor <- m[m$sample!='Sham4',];Xdrop <- model.matrix(~batch+condition,without_anchor)
audit <- data.frame(design=c('all_seven_samples','omit_only_batch4_sham'),n_samples=c(nrow(X),nrow(Xdrop)),n_columns=c(ncol(X),ncol(Xdrop)),rank=c(qr(X)$rank,qr(Xdrop)$rank))
write.table(audit,file.path(out,'fusion_design_rank.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Median-of-ratios size factors across the seven samples, fixed for descriptive deletion checks.
geomean <- exp(rowMeans(log(counts)));valid <- is.finite(geomean)&geomean>0
sf <- apply(counts[valid,,drop=FALSE]/geomean[valid],2,median); sf <- sf/exp(mean(log(sf)))
norm <- sweep(counts,2,sf,'/');logn <- log2(norm+1)
ko <- m$condition=='MYMK_KO'; sham4 <- m$sample=='Sham4'
within <- rowMeans(logn[,ko,drop=FALSE])-logn[,sham4]
full <- read.delim('results/lineage_gene_statistics.tsv')
full <- full[full$contrast=='overload_fusion_blockade_vs_sham_batch_adjusted',]
a <- match(full$gene,rownames(counts));stopifnot(!anyNA(a))
res <- data.frame(gene=full$gene,original_lfc=full$log2FoldChange,original_padj=full$padj,within_batch_descriptive_log2_difference=within[a])
for(i in which(ko)) res[[paste0('omit_',m$sample[i])]] <- rowMeans(logn[a,setdiff(which(ko),i),drop=FALSE])-logn[a,sham4]
res$batch4_sham_minus_batch3_sham <- logn[a,sham4]-rowMeans(logn[a,m$condition=='SHAM'&m$batch=='3',drop=FALSE])
write.table(res,file.path(out,'fusion_within_batch_gene_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
sig <- !is.na(res$original_padj)&res$original_padj<.05
stats <- data.frame(n_original_significant=sum(sig),n_direction_concordant=sum(sign(res$original_lfc[sig])==sign(res$within_batch_descriptive_log2_difference[sig])))
# Explicit rowwise check, avoiding any implicit recycling of effect directions.
stats$n_direction_retained_all_KO_deletions <- sum(rowSums(sign(as.matrix(res[sig,grep('^omit_',names(res)),drop=FALSE]))==sign(res$original_lfc[sig]))==sum(grepl('^omit_',names(res))))
stats$spearman_all_tested <- cor(res$original_lfc,res$within_batch_descriptive_log2_difference,method='spearman',use='complete.obs')
write.table(stats,file.path(out,'fusion_sensitivity_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(res[res$gene %in% c('Mymk','Runx1'),],file.path(out,'fusion_selected_gene_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file='results/R_session_fusion.txt')
print(audit);print(stats);print(res[res$gene %in% c('Mymk','Runx1'),])
