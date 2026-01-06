rm(list = ls())

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/wu-gu-prostate-bladder/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(magrittr)
library(patchwork)
library(lme4)
library(presto)

# Data 
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta_amended.csv", row.names = 1)
norm <- norm[rownames(norm) != "NegProbe-WTX", rownames(meta)]
cts <- cts[rownames(cts) != "NegProbe-WTX", rownames(meta)]

# This analysis only looks at the prostate data
prostate <- meta[(meta$cancer_type == "prostate"),] |> rownames()
cts <- cts[,prostate]
norm <- norm[,prostate]
meta <- meta[prostate,]

# For visualization
snorm <- apply(X = norm, MARGIN = 1, FUN = scale) |> t()
colnames(snorm) <- sprintf("%03d", meta$roi)
bc <- limma::removeBatchEffect(x = norm, batch = meta$patient_deid, design = model.matrix(~sub_types_v2, data = meta))
sbc <- apply(X = bc, MARGIN = 1, FUN = scale) |> t()
colnames(sbc) <- sprintf("%03d", meta$roi)

# From 7/18/25 email -----------------------------------------------------------

## Request 3: Page 9 of the attached PDF is from the preliminary analysis during
## our first meeting. Can we get the PDF files with diagram and heatmap 
## (including batch-corrected) for this comparison? Thank you!
idx <- meta$sub_types %in% c("Ductal", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- (mm[dds@colData$sub_types == "Ductal", ] |> colMeans()) - 
  (mm[dds@colData$sub_types == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 105
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.8499692

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "Ductal", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "03 - DESeq2_Ductal_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "03 - DESeq2_Ductal_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Ductal"="yellow2", "Acinar_crib"="green4")))
pdf(file = "03 - DESeq2_Ductal_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 8)
Heatmap(matrix = mat, top_annotation = ha, 
        column_split = meta[idx,]$sub_types_v2,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "03 - DESeq2_Ductal_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 8)
Heatmap(matrix = mat, top_annotation = ha, 
        column_split = meta[idx,]$sub_types_v2,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


## Request 10: Is there a batch-corrected version of the file below?
### 10.4 - DESeq2_Ductal_vs_Acinar_core_paired_top_genes_heatmap

## This got tricky, since we had to pseudobulk to keep the paired structure.
meta[meta$patient_deid == "pt_30", "sub_types_v2"]

# [1] "Acinar_non_crib_LG" "Acinar_crib"        "Acinar IDC-P_crib"  "Ductal"            

## See, there are two AOIs from pt_30 that are in the Acinar group of interest
## that can be paired with a Ductal AOI from the same patient.

## Ting wants the batch-corrected counts for visualization, which essentially
## subtracts out the patient differences.

## So, I have to bulk these two AOIs together, then do batch correction, then
## visualize again.

meta$sub_types_v4 <- ifelse(test = meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG", "Acinar_crib"), 
                            yes = "Acinar", no = meta$sub_types_v2)

pairs <- table(meta$tma_core_number, meta$sub_types_v4) |> 
  as.data.frame.matrix() |> 
  select(Ductal, Acinar)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

which((pairs > 1) |> rowSums() > 0) |> names() # This gets the core for which we will bulk
# [1] "TMA1D1"

meta_copy <- meta
meta_copy$roi %<>% sprintf("%03d", .)
meta_copy[meta_copy$tma_core_number == "TMA1D1" & meta_copy$sub_types_v4 == "Acinar",]$roi <- 
  meta_copy[meta_copy$tma_core_number == "TMA1D1" & meta_copy$sub_types_v4 == "Acinar",]$roi |> as.list() |> do.call(what = paste, args = _)
meta_copy$roi %<>% gsub(pattern = " ", replacement = "_", x = .)

library(data.table)
psb <- presto::collapse_counts(counts_mat = cts |> as.matrix(), # pseudo-bulking
                               meta_data = meta_copy, 
                               varnames = c("roi", "patient_deid", "tma_core_number", "sub_types_v4"))
qs <- apply(X = psb$counts_mat, MARGIN = 2, FUN = quantile, 0.75)
psb$meta_data$q_norm_qfactors <- (qs / EnvStats::geoMean(qs))
colnames(psb$counts_mat) <- psb$meta_data$roi
rownames(psb$meta_data) <- psb$meta_data$roi

norm_tmp <- log2((psb$counts_mat / qs) + 1)
snorm_tmp <- apply(X = norm_tmp, MARGIN = 1, FUN = scale) |> t()
colnames(snorm_tmp) <- colnames(norm_tmp)
bc_tmp <- limma::removeBatchEffect(x = norm_tmp, batch = psb$meta_data$patient_deid, design = model.matrix(~sub_types_v4, data = psb$meta_data))
sbc_tmp <- apply(X = bc_tmp, MARGIN = 1, FUN = scale) |> t()
colnames(sbc_tmp) <- colnames(bc_tmp)

paired_meta <- psb$meta_data[psb$meta_data$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v4 %in% c("Ductal", "Acinar"))
paired_meta$group <- paired_meta$sub_types_v4

dds <- DESeq2::DESeqDataSetFromMatrix(countData = psb$counts_mat[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Acinar", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 138
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.86309

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "Ductal", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar", hjust = 1, fontface = "bold")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
coreidx <- paired_meta |> arrange(group, tma_core_number) |> rownames()
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
ha <- HeatmapAnnotation(sub_type = paired_meta[coreidx,]$group,
                        tma_core = paired_meta[coreidx,]$tma_core_number,
                        col = list("sub_type"=c("Acinar"="orange", "Ductal"="yellow2"), 
                                   "tma_core"=corecols))

mat <- snorm_tmp[toplot, coreidx]
Heatmap(matrix = mat, column_split = paired_meta[coreidx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)

mat <- sbc_tmp[toplot, coreidx]
Heatmap(matrix = mat, column_split = paired_meta[coreidx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)

# It looks terrible. I do not think it really makes sense. I am going to 
# suggest against doing this. I tried.


## Request 11: Re-make the heatmap, based on effect size only. Add the 
## annotation she requested.
meta$sub_types_v5 <- case_when(meta$roi %in% c(48, 49) ~ "Precursor", 
                               meta$sub_types_v2 %in% c("Acinar IDC-P", "Acinar IDC-P_crib", "AIP") ~ "Intraductal_spread", 
                               T ~ meta$sub_types_v2)

small_meta <- meta |> filter(sub_types_v5 %in% c("Precursor", "Intraductal_spread"))
small_meta$group <- small_meta$sub_types_v5

# Calculating the difference in means across the two groups for every gene.
mm <- model.matrix(~0+group, data = small_meta)
small_norm <- norm[,rownames(small_meta)]
groupmeans <- (as.matrix(small_norm) %*% mm) |> sweep(x = _, MARGIN = 2, STATS = colSums(mm), FUN = "/")
colnames(groupmeans) <- gsub(pattern = "group", replacement = "", x = colnames(groupmeans))
groupmeans <- as.data.frame(groupmeans)

toplot <- small_norm[which((abs(groupmeans[,"Precursor"]-groupmeans[,"Intraductal_spread"])) > 2),] |> rownames()
idx <- sprintf("%03d", small_meta$roi)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(group = small_meta$group,
                        sub_type = small_meta$sub_types_v2,
                        col = list("group"=c("Precursor"="dodgerblue", "Intraductal_spread"="gold3"), 
                                   "sub_type"=c("Acinar IDC-P"="cyan", "Acinar IDC-P_crib"="blue", "AIP"="purple")))
pdf(file = "11 - DESeq2_Precursor_vs_Intraductal_spread_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = small_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "11 - DESeq2_Precursor_vs_Intraductal_spread_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = small_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

# The batch-corrected does not look great. Also, I think it makes sense to do a
# ssGSEA thing. Ting will need to choose the pathways.

## Request 12: Sorry Cole, we decide not to use 2 AOIs named as “Acinar IDC-P” 
## (AOI 14 and 18; these are non-cribriform Acinar IDC-P, although I didn’t 
## indicate in the name), because the morphology in these 2 AOIs is a bit 
## contraversial. Dense cribriform is more universally accepted as IDC-P. 
## Therefore, we decide to use “Acinar IDC-P_crib” to represent all Acinar IDC-P.
## I wonder if you could re-run 14.1, 14.2, 14.3, 14.4 and 16 comparisons in the 
## folders below, after deleting “Acinar IDC-P” (AOI 14 and 18), and use 
## “Acinar IDC-P_crib” to compare with other groups. Thank you very much!

## 14.1) All Acinar IDC-P vs all acinar carcinoma ------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar_non_crib_HG", "Acinar_crib") ~ "all_acinar_carcinoma",
  T ~ meta$sub_types_v2
)

idx <- ((meta$group %in% c("all_acinar_IDC-P", "all_acinar_carcinoma")) & (!meta$roi %in% c(14, 18)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta[idx,]$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "all_acinar_carcinoma", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 97
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.043615

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "all_acinar_carcinoma", hjust = 1, fontface = "bold")
ggsave(filename = "14.1 - DESeq2_all_acinar_IDC-P_vs_all_acinar_carcinoma.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "14.1 - DESeq2_all_acinar_IDC-P_vs_all_acinar_carcinoma.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("Acinar IDC-P_crib"="lightblue", "Acinar IDC-P"="dodgerblue", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "14.1 - DESeq2_all_acinar_IDC-P_vs_all_acinar_carcinoma_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "14.1 - DESeq2_all_acinar_IDC-P_vs_all_acinar_carcinoma_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 14.2) All Acinar IDC-P vs Acinar non-crib_LG --------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  T ~ meta$sub_types_v2
)
idx <- (meta$group %in% c("all_acinar_IDC-P", "Acinar_non_crib_LG") & (!meta$roi %in% c(14, 18)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib_LG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 110
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9661684

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_LG", hjust = 1, fontface = "bold")
ggsave(filename = "14.2 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_LG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "14.2 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_LG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("Acinar IDC-P_crib"="lightblue", "Acinar IDC-P"="dodgerblue", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "14.2 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_LG_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "14.2 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_LG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 14.3) All Acinar IDC-P vs Acinar non-crib_HG --------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  T ~ meta$sub_types_v2
)
idx <- (meta$group %in% c("all_acinar_IDC-P", "Acinar_non_crib_HG") & (!meta$roi %in% c(14, 18)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib_HG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 227
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.553474

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_HG", hjust = 1, fontface = "bold")
ggsave(filename = "14.3 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_HG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "14.3 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_HG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("Acinar IDC-P_crib"="lightblue", "Acinar IDC-P"="dodgerblue", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "14.3 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_HG_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "14.3 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_HG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 14.4) All Acinar IDC-P vs Acinar crib ---------------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  T ~ meta$sub_types_v2
)
idx <- (meta$group %in% c("all_acinar_IDC-P", "Acinar_crib") & (!meta$roi %in% c(14, 18)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 6
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.149369

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "14.4 - DESeq2_all_acinar_IDC-P_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "14.4 - DESeq2_all_acinar_IDC-P_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("Acinar IDC-P_crib"="lightblue", "Acinar IDC-P"="dodgerblue", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "14.4 - DESeq2_all_acinar_IDC-P_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "14.4 - DESeq2_all_acinar_IDC-P_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 16) All Acinar IDC-P vs AIP -------------------------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  T ~ meta$sub_types_v2
)
idx <- (meta$group %in% c("all_acinar_IDC-P", "AIP") & (!meta$roi %in% c(14, 18)))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "AIP", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 55
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.7380318

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "AIP", hjust = 1, fontface = "bold")
ggsave(filename = "16 - DESeq2_all_acinar_IDC-P_vs_AIP.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "16 - DESeq2_all_acinar_IDC-P_vs_AIP.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("AIP"="purple", "Acinar IDC-P"="cyan", "Acinar IDC-P_crib"="blue")))
pdf(file = "16 - DESeq2_all_acinar_IDC-P_vs_AIP_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,] == "AIP",
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "16 - DESeq2_all_acinar_IDC-P_vs_AIP_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,] == "AIP",
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


## Request 13: Could you please add 4 paired comparisons (within the same tumor) 
## as listed below? Thank you! I’ve put the corresponding AOIs and names in the 
## excel. Please let me know if you have questions.
### Paired acinar IDC-P_crib vs. all Acinar Ca
### Paired acinar IDC-P_crib vs. Acinar_non-crib_LG
### Paired acinar IDC-P_crib vs. Acinar_non-crib_HG
### Paired acinar IDC-P_crib vs. Acinar_crib

## 13.1) Paired acinar IDC-P_crib vs all Acinar
meta$sub_types_v6 <- ifelse(test = meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG", "Acinar_crib"), 
                            yes = "all_acinar_carcinoma", no = meta$sub_types_v2)

pairs <- table(meta$tma_core_number, meta$sub_types_v6) |> 
  as.data.frame.matrix() |> 
  select(all_acinar_carcinoma, `Acinar IDC-P_crib`)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v6 %in% c("all_acinar_carcinoma", "Acinar IDC-P_crib"))
paired_meta$group <- paired_meta$sub_types_v6

library(data.table)
psb <- presto::collapse_counts(counts_mat = cts[,rownames(paired_meta)] |> as.matrix(), # pseudo-bulking
                               meta_data = paired_meta, 
                               varnames = c("group", "tma_core_number"))
qs <- apply(X = psb$counts_mat, MARGIN = 2, FUN = quantile, 0.75) # Normalizing the pseudo-bulked data
psb$meta_data$q_norm_qfactors <- (qs / EnvStats::geoMean(qs))
psb$meta_data$patient_deid <- plyr::mapvalues(x = psb$meta_data$tma_core_number, from = meta$tma_core_number, to = meta$patient_deid)
roimap <- paired_meta |> 
  dplyr::group_by(group, tma_core_number) |> 
  dplyr::mutate(roi = sprintf("%03d", roi)) |>
  dplyr::reframe(rois = stringr::str_c(roi, collapse = "_"))
roimap$id <- paste(roimap$group, roimap$tma_core_number, sep = "_")
psb$meta_data$id <- paste(psb$meta_data$group, psb$meta_data$tma_core_number, sep = "_")
psb$meta_data$roi <- plyr::mapvalues(x = psb$meta_data$id, from = roimap$id, to = roimap$rois)
colnames(psb$counts_mat) <- psb$meta_data$roi
rownames(psb$meta_data) <- psb$meta_data$roi

dds <- DESeq2::DESeqDataSetFromMatrix(countData = psb$counts_mat, colData = psb$meta_data, design = ~group+tma_core_number)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_carcinoma", ] |> colMeans()) - (mm[dds@colData$group == "Acinar IDC-P_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 179
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.7911838

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "all_acinar_carcinoma", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar IDC-P_crib", hjust = 1, fontface = "bold")
ggsave(filename = "13 - DESeq2_all_acinar_carcinoma_vs_Acinar_IDC-P_crib_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "13 - DESeq2_all_acinar_carcinoma_vs_Acinar_IDC-P_crib_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(psb$meta_data$tma_core_number), palette = "brewers")
coreidx <- psb$meta_data |> arrange(group, tma_core_number) |> rownames()
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- (log2(DESeq2::counts(dds, normalized = T)+1)[toplot, coreidx] |> apply(MARGIN = 1, FUN = scale) |> t())
colnames(mat) <- psb$meta_data[coreidx,]$roi
ha <- HeatmapAnnotation(sub_type = psb$meta_data[coreidx,]$group,
                        tma_core = psb$meta_data[coreidx,]$tma_core_number,
                        col = list("sub_type"=c("all_acinar_carcinoma"="orange", "Acinar IDC-P_crib"="blue"), 
                                   "tma_core"=corecols))
pdf(file = "13 - DESeq2_all_acinar_carcinoma_vs_Acinar_IDC-P_crib_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = psb$meta_data[coreidx,]$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


## 13.2) Paired acinar IDC-P_crib vs Acinar_non-crib_LG
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Acinar_non_crib_LG, `Acinar IDC-P_crib`)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar IDC-P_crib"))
paired_meta$group <- paired_meta$sub_types_v2

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Acinar_non_crib_LG", ] |> colMeans()) - (mm[dds@colData$group == "Acinar IDC-P_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 74
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9538702

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "Acinar_non_crib_LG", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar IDC-P_crib", hjust = 1, fontface = "bold")
ggsave(filename = "13 - DESeq2_Acinar_non_crib_LG_vs_Acinar_IDC-P_crib_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "13 - DESeq2_Acinar_non_crib_LG_vs_Acinar_IDC-P_crib_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- paired_meta$roi |> sprintf("%03d", ... = _)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = paired_meta$group,
                        tma_core = paired_meta$tma_core_number,
                        col = list("sub_type"=c("Acinar_non_crib_LG"="pink", "Acinar IDC-P_crib"="blue"), 
                                   "tma_core"=corecols))
pdf(file = "13 - DESeq2_Acinar_non_crib_LG_vs_Acinar_IDC-P_crib_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "13 - DESeq2_Acinar_non_crib_LG_vs_Acinar_IDC-P_crib_core_paired_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


## 13.3) Paired acinar IDC-P_crib vs Acinar_non-crib_HG
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Acinar_non_crib_HG, `Acinar IDC-P_crib`)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar IDC-P_crib"))
paired_meta$group <- paired_meta$sub_types_v2

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Acinar_non_crib_HG", ] |> colMeans()) - (mm[dds@colData$group == "Acinar IDC-P_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 163
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.8343528

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "Acinar_non_crib_HG", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar IDC-P_crib", hjust = 1, fontface = "bold")
ggsave(filename = "13 - DESeq2_Acinar_non_crib_HG_vs_Acinar_IDC-P_crib_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "13 - DESeq2_Acinar_non_crib_HG_vs_Acinar_IDC-P_crib_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- paired_meta$roi |> sprintf("%03d", ... = _)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = paired_meta$group,
                        tma_core = paired_meta$tma_core_number,
                        col = list("sub_type"=c("Acinar_non_crib_HG"="firebrick", "Acinar IDC-P_crib"="blue"), 
                                   "tma_core"=corecols))
pdf(file = "13 - DESeq2_Acinar_non_crib_HG_vs_Acinar_IDC-P_crib_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "13 - DESeq2_Acinar_non_crib_HG_vs_Acinar_IDC-P_crib_core_paired_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 13.4) Paired acinar IDC-P_crib vs Acinar_crib
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Acinar_crib, `Acinar IDC-P_crib`)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Acinar_crib", "Acinar IDC-P_crib"))
paired_meta$group <- paired_meta$sub_types_v2

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Acinar_crib", ] |> colMeans()) - (mm[dds@colData$group == "Acinar IDC-P_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 0
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9409709

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "Acinar_crib", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar IDC-P_crib", hjust = 1, fontface = "bold")
ggsave(filename = "13 - DESeq2_Acinar_crib_vs_Acinar_IDC-P_crib_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "13 - DESeq2_Acinar_crib_vs_Acinar_IDC-P_crib_core_paired.csv")


# From 7/23/25 email -----------------------------------------------------------

# Data 
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta_amended.csv", row.names = 1)
norm <- norm[rownames(norm) != "NegProbe-WTX", rownames(meta)]
cts <- cts[rownames(cts) != "NegProbe-WTX", rownames(meta)]

# This analysis only looks at the prostate data
bladder <- meta[(meta$cancer_type == "bladder"),] |> rownames()
cts <- cts[,bladder]
norm <- norm[,bladder]
meta <- meta[bladder,]

# For visualization
snorm <- apply(X = norm, MARGIN = 1, FUN = scale) |> t()
colnames(snorm) <- sprintf("%03d", meta$roi)
bc <- limma::removeBatchEffect(x = norm, batch = meta$patient_deid, design = model.matrix(~sub_types_v2, data = meta))
sbc <- apply(X = bc, MARGIN = 1, FUN = scale) |> t()
colnames(sbc) <- sprintf("%03d", meta$roi)

## Request 2: Please delete AOI 31 (SM_solid) for the comparisons below. 
## Request 3: Could you please compare the following like what we did for 
## prostate project?

### All MP vs. all CIUC
idx <- ((meta$sub_types_v2 %in% c("MP", "CIUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "MP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "CIUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 268
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.434091

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "MP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "CIUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_MP_vs_CIUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_MP_vs_CIUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("MP"="gold", "CIUC"="hotpink")))
pdf(file = "bladder - DESeq2_MP_vs_CIUC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_MP_vs_CIUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All PC vs. all CIUC
idx <- ((meta$sub_types_v2 %in% c("PC", "CIUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "PC", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "CIUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 106
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
#[1] 1.167703

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "PC", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "CIUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_PC_vs_CIUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_PC_vs_CIUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("PC"="grey", "CIUC"="hotpink")))
pdf(file = "bladder - DESeq2_PC_vs_CIUC_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_PC_vs_CIUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All SM (exclude AOI 31) vs. all CIUC
idx <- ((meta$sub_types_v2 %in% c("SM", "CIUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "SM", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "CIUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 596
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.010728

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "SM", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "CIUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_SM_vs_CIUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_SM_vs_CIUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("SM"="grey2", "CIUC"="hotpink")))
pdf(file = "bladder - DESeq2_SM_vs_CIUC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_SM_vs_CIUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All MP vs. all PUC
idx <- ((meta$sub_types_v2 %in% c("MP", "PUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "MP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "PUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 357
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.8531702

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "MP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "PUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_MP_vs_PUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_MP_vs_PUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("MP"="gold", "PUC"="limegreen")))
pdf(file = "bladder - DESeq2_MP_vs_PUC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_MP_vs_PUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All PC vs. all PUC
idx <- ((meta$sub_types_v2 %in% c("PC", "PUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "PC", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "PUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 290
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.7091075

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "PC", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "PUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_PC_vs_PUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_PC_vs_PUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("PC"="grey", "PUC"="limegreen")))
pdf(file = "bladder - DESeq2_PC_vs_PUC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_PC_vs_PUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All SM (exclude AOI 31) vs. all PUC
idx <- ((meta$sub_types_v2 %in% c("SM", "PUC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "SM", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "PUC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 315
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.6206226

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "SM", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "PUC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_SM_vs_PUC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_SM_vs_PUC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("SM"="grey2", "PUC"="limegreen")))
pdf(file = "bladder - DESeq2_SM_vs_PUC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_SM_vs_PUC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

### All MP vs. All PC
idx <- ((meta$sub_types_v2 %in% c("MP", "PC")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "MP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "PC", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 30
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.180628

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "MP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "PC", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_MP_vs_PC.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_MP_vs_PC.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("MP"="gold", "PC"="grey")))
pdf(file = "bladder - DESeq2_MP_vs_PC_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_MP_vs_PC_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All MP vs. all SM (exclude AOI 31)
idx <- ((meta$sub_types_v2 %in% c("MP", "SM")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "MP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "SM", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 925
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.554779

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "MP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "SM", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_MP_vs_SM.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_MP_vs_SM.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("MP"="gold", "SM"="grey2")))
pdf(file = "bladder - DESeq2_MP_vs_SM_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_MP_vs_SM_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()


### All PC vs. all SM (exclude AOI 31)
idx <- ((meta$sub_types_v2 %in% c("PC", "SM")) & (meta$roi != 31))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2)
dds@colData$sizeFactor <- dds$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "PC", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "SM", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 543
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.374347

res$target <- rownames(res)
ggplot() + 
  scattermore::geom_scattermore(data = res, mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 4, color = "grey") + 
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) < 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "orange") +
  scattermore::geom_scattermore(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], mapping = aes(x = log2FoldChange, y = -log10(padj)), pointsize = 3, color = "red") + 
  ggthemes::theme_par() + 
  geom_hline(yintercept = -log10(0.05), color = "black", linewidth = 0.25, linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1), color = "black", linewidth = 0.25, linetype = "dashed") + 
  ggrepel::geom_text_repel(data = res[res$padj < 0.05 & abs(res$log2FoldChange) > 1,], 
                           mapping = aes(x = log2FoldChange, y = -log10(padj), label = target), 
                           color = "black", min.segment.length = 0, box.padding = 0.25, max.overlaps = 30, size = 3, segment.size = 0.25) + 
  annotate(geom = "segment", x = 0, xend = 0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) + 
  annotate(geom = "text", x = 0.6, y = -2, label = "PC", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "SM", hjust = 1, fontface = "bold")
ggsave(filename = "bladder - DESeq2_PC_vs_SM.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "bladder - DESeq2_PC_vs_SM.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("PC"="grey", "SM"="grey2")))
pdf(file = "bladder - DESeq2_PC_vs_SM_top_genes_heatmap.pdf", width = 8, height = 18)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

pdf(file = "bladder - DESeq2_PC_vs_SM_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 18)
mat <- sbc[toplot, idx]
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
        row_names_gp = gpar(fontsize = 3), 
        show_column_names = T, column_names_gp = gpar(fontsize = 10),
        cluster_columns = F,
        name = "Batch-corrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

### Session --------------------------------------------------------------------
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] presto_1.0.0          data.table_1.16.4     Rcpp_1.0.14           lme4_1.1-37           Matrix_1.7-3          patchwork_1.3.0       magrittr_2.0.3       
# [8] ComplexHeatmap_2.22.0 ggplot2_3.5.2         dplyr_1.1.4          
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4                deldir_2.0-4                gridExtra_2.3               rlang_1.1.5                 clue_0.3-66                
# [6] GetoptLong_1.0.5            matrixStats_1.5.0           compiler_4.4.2              spatstat.geom_3.3-6         png_0.1-8                  
# [11] vctrs_0.6.5                 stringr_1.5.1               pkgconfig_2.0.3             shape_1.4.6.1               crayon_1.5.3               
# [16] XVector_0.46.0              labeling_0.4.3              nloptr_2.2.1                UCSC.utils_1.2.0            purrr_1.0.4                
# [21] zlibbioc_1.52.0             GenomeInfoDb_1.42.3         jsonlite_2.0.0              EnvStats_3.0.0              SnowballC_0.7.1            
# [26] DelayedArray_0.32.0         spatstat.utils_3.1-3        BiocParallel_1.40.2         irlba_2.3.5.1               parallel_4.4.2             
# [31] cluster_2.1.6               R6_2.6.1                    stringi_1.8.4               RColorBrewer_1.1-3          spatstat.data_3.1-6        
# [36] limma_3.62.2                reticulate_1.42.0           spatstat.univar_3.1-2       boot_1.3-31                 scattermore_1.2            
# [41] GenomicRanges_1.58.0        SummarizedExperiment_1.36.0 iterators_1.0.14            IRanges_2.40.1              splines_4.4.2              
# [46] tidyselect_1.2.1            viridis_0.6.5               rstudioapi_0.17.1           abind_1.4-8                 doParallel_1.0.17          
# [51] codetools_0.2-20            plyr_1.8.9                  lattice_0.22-6              tibble_3.2.1                Biobase_2.66.0             
# [56] withr_3.0.2                 askpass_1.2.1               polyclip_1.10-7             circlize_0.4.16             mclust_6.1.1               
# [61] pillar_1.10.1               lsa_0.73.3                  BiocManager_1.30.25         MatrixGenerics_1.18.1       renv_1.1.1                 
# [66] foreach_1.5.2               stats4_4.4.2                reformulas_0.4.0            generics_0.1.3              S4Vectors_0.44.0           
# [71] munsell_0.5.1               scales_1.3.0                minqa_1.2.8                 glue_1.8.0                  InSituType_2.0             
# [76] tools_4.4.2                 RSpectra_0.16-2             locfit_1.5-9.12             rbibutils_2.3               umap_0.2.10.0              
# [81] colorspace_2.1-1            SingleCellExperiment_1.28.1 nlme_3.1-166                GenomeInfoDbData_1.2.13     cli_3.6.4                  
# [86] ggthemes_5.1.0              S4Arrays_1.6.0              viridisLite_0.4.2           uwot_0.2.3                  gtable_0.3.6               
# [91] DESeq2_1.46.0               digest_0.6.37               BiocGenerics_0.52.0         SparseArray_1.6.2           ggrepel_0.9.6              
# [96] rjson_0.2.23                farver_2.1.2                lifecycle_1.0.4             httr_1.4.7                  GlobalOptions_0.1.2        
# [101] statmod_1.5.0               openssl_2.3.2               MASS_7.3-61                

