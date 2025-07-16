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

## 7) Acinar non-crib vs Acinar crib -------------------------------------------
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar_non_crib_HG", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  ((mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans()) + 
  (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans()))/2 - 
  (mm[dds@colData$sub_types_v2 == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 149
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9432348

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
  annotate(geom = "text", x = 0.6, y = -2, label = "Acinar_non_crib", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "07 - DESeq2_Acinar_non_crib_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "07 - DESeq2_Acinar_non_crib_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Acinar_crib"="green4")))
pdf(file = "07 - DESeq2_Acinar_non_crib_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG"),
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
pdf(file = "07 - DESeq2_Acinar_non_crib_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG"),
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

## 7) Acinar non-crib HG vs Acinar crib ----------------------------------------
contr <- 
  (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans()) - 
  (mm[dds@colData$sub_types_v2 == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 143
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# 0.718904

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "07 - DESeq2_Acinar_non_crib_HG_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "07 - DESeq2_Acinar_non_crib_HG_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_crib")
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Acinar_crib"="green4")))
pdf(file = "07 - DESeq2_Acinar_non_crib_HG_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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
pdf(file = "07 - DESeq2_Acinar_non_crib_HG_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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

## 7) Acinar non-crib LG vs Acinar crib ----------------------------------------
contr <- 
  (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans()) - 
  (mm[dds@colData$sub_types_v2 == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 105
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# 1.214444

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "07 - DESeq2_Acinar_non_crib_LG_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "07 - DESeq2_Acinar_non_crib_LG_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar_crib")
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Acinar_crib"="green4")))
pdf(file = "07 - DESeq2_Acinar_non_crib_LG_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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
pdf(file = "07 - DESeq2_Acinar_non_crib_LG_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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

## 8) Ductal vs Acinar non-crib  -----------------------------------------------
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar_non_crib_HG", "Ductal")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "Ductal", ] |> colMeans()) -
  ((mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans()) + 
     (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans()))/2
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 511
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9110375

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib", hjust = 1, fontface = "bold")
ggsave(filename = "08 - DESeq2_Ductal_vs_Acinar_non_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "08 - DESeq2_Ductal_vs_Acinar_non_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Ductal"="yellow2")))
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG"),
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
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG"),
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

## 8) Ductal vs Acinar non-crib HG ---------------------------------------------
contr <- 
  (mm[dds@colData$sub_types_v2 == "Ductal", ] |> colMeans()) -
  (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans()) 
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 490
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.210526

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_HG", hjust = 1, fontface = "bold")
ggsave(filename = "08 - DESeq2_Ductal_vs_Acinar_non_crib_HG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_HG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Ductal")
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Ductal"="yellow2")))
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_HG_top_genes_heatmap.pdf", width = 8, height = 16)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_HG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 16)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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

## 8) Ductal vs Acinar non-crib LG ---------------------------------------------
contr <- 
  (mm[dds@colData$sub_types_v2 == "Ductal", ] |> colMeans()) -
  (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans()) 
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 270
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.7356175

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_LG", hjust = 1, fontface = "bold")
ggsave(filename = "08 - DESeq2_Ductal_vs_Acinar_non_crib_LG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_LG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Ductal")
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2, 
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Ductal"="yellow2")))
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_LG_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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
pdf(file = "08 - DESeq2_Ductal_vs_Acinar_non_crib_LG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = meta[idx,]$sub_types_v2,
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

## 9) Ductal vs. Acinar crib ---------------------------------------------------
# Already covered in the original DE analysis.

## 10.1) Paired Ductal vs. Acinar non-crib -------------------------------------
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Ductal, Acinar_non_crib_LG, Acinar_non_crib_HG)
pairs <- pairs[rowSums(pairs) > 1,]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Ductal", "Acinar_non_crib_LG", "Acinar_non_crib_HG"))
paired_meta$roi <- sprintf("%03d", paired_meta$roi)
paired_meta$group <- paired_meta$sub_types_v2 |> gsub(pattern = "_[A-Z]{2}", replacement = "")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- paired_meta$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 117
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# 1.724201

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib", hjust = 1, fontface = "bold")
ggsave(filename = "10.1 - DESeq2_Ductal_vs_Acinar_non_crib_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "10.1 - DESeq2_Ductal_vs_Acinar_non_crib_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
coreidx <- paired_meta |> arrange(group, tma_core_number) |> pull(roi)
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, coreidx]
paired_meta <- paired_meta |> arrange(group, tma_core_number)
ha <- HeatmapAnnotation(sub_type = paired_meta$sub_types_v2,
                        tma_core = paired_meta$tma_core_number,
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Ductal"="yellow2"), 
                                   "tma_core"=corecols))
pdf(file = "10.1 - DESeq2_Ductal_vs_Acinar_non_crib_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, coreidx]
pdf(file = "10.1 - DESeq2_Ductal_vs_Acinar_non_crib_core_paired_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = paired_meta$group,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 10.2) Paired Ductal vs. Acinar non-crib LG ----------------------------------
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Ductal, Acinar_non_crib_LG)
pairs <- pairs[rowSums(pairs) > 1,]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Ductal", "Acinar_non_crib_LG"))
paired_meta$group <- paired_meta$sub_types_v2

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- paired_meta$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib_LG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 17
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.55783

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_LG", hjust = 1, fontface = "bold")
ggsave(filename = "10.2 - DESeq2_Ductal_vs_Acinar_non_crib_LG_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "10.2 - DESeq2_Ductal_vs_Acinar_non_crib_LG_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
coreidx <- paired_meta |> arrange(group, tma_core_number) |> pull(roi)
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, coreidx]
paired_meta <- paired_meta |> arrange(group, tma_core_number)
ha <- HeatmapAnnotation(sub_type = paired_meta$sub_types_v2,
                        tma_core = paired_meta$tma_core_number,
                        col = list("sub_type"=c("Acinar_non_crib_HG"="red3", "Acinar_non_crib_LG"="pink", "Ductal"="yellow2"), 
                                   "tma_core"=corecols))
pdf(file = "10.2 - DESeq2_Ductal_vs_Acinar_non_crib_LG_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = paired_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, coreidx]
pdf(file = "10.2 - DESeq2_Ductal_vs_Acinar_non_crib_LG_core_paired_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, top_annotation = ha, column_split = paired_meta$group,
        row_names_gp = gpar(fontsize = 7), 
        column_names_gp = gpar(fontsize = 7),
        cluster_columns = F,
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

## 10.3) Paired Ductal vs. Acinar HG -------------------------------------------
meta$sub_types_v3 <- ifelse(test = meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_crib"), 
                                     yes = "Acinar_HG", no = meta$sub_types_v2)
pairs <- table(meta$tma_core_number, meta$sub_types_v3) |> 
  as.data.frame.matrix() |> 
  select(Ductal, Acinar_HG)
pairs <- pairs[apply(X = pairs, MARGIN = 1, FUN = GeomxTools::ngeoMean) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v3 %in% c("Ductal", "Acinar_HG"))
paired_meta$group <- paired_meta$sub_types_v3

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- paired_meta$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_HG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 1
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.767442

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_HG", hjust = 1, fontface = "bold")
ggsave(filename = "10.3 - DESeq2_Ductal_vs_Acinar_HG_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "10.3 - DESeq2_Ductal_vs_Acinar_HG_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(paired_meta$tma_core_number), palette = "brewers")
coreidx <- paired_meta |> arrange(group, tma_core_number) |> rownames()
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)

paired_meta <- paired_meta |> arrange(group, tma_core_number)
paired_meta$MGP <- norm["MGP", coreidx] |> t()
paired_meta[coreidx,]$roi <- sprintf("%03d", paired_meta$roi)
ggplot(data = paired_meta) + 
  stat_summary(mapping = aes(x = group, y = MGP, color = group), geom = "crossbar", fun = mean) + 
  scale_color_manual(values = c("Acinar_HG"="violet", "Ductal"="yellow2")) +
  ggnewscale::new_scale_color() +
  geom_line(mapping = aes(x = group, y = MGP, group = tma_core_number), color = "grey") +
  geom_point(mapping = aes(x = group, y = MGP, color = tma_core_number), size = 4, shape = 21, stroke = 2, fill = "white") + 
  scale_color_manual(values = corecols) +
  ggrepel::geom_text_repel(mapping = aes(x = group, y = MGP, label = roi), point.padding = 10, size = 3) +
  ggthemes::theme_par() + 
  theme(axis.title.x = element_blank()) +
  labs(title = "MGP", y = "Normalized Expression", color = "TMA core")
ggsave(filename = "10.3 - DESeq2_Ductal_vs_Acinar_HG_core_paired_MGP_expression.pdf", height = 6, width = 6)

## 10.4) Paired Ductal vs. Acinar ----------------------------------------------

# * Note that this one is especially difficult because we need to pseudo-bulk *

meta$sub_types_v4 <- ifelse(test = meta$sub_types_v2 %in% c("Acinar_non_crib_HG", "Acinar_non_crib_LG", "Acinar_crib"), 
                                     yes = "Acinar", no = meta$sub_types_v2)
pairs <- table(meta$tma_core_number, meta$sub_types_v4) |> 
  as.data.frame.matrix() |> 
  select(Ductal, Acinar)
pairs <- pairs[rowMeans(pairs > 0) == 1, ]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v4 %in% c("Ductal", "Acinar"))
paired_meta$group <- paired_meta$sub_types_v4

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
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Acinar", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 262
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.880692

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
ggsave(filename = "10.4 - DESeq2_Ductal_vs_Acinar_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "10.4 - DESeq2_Ductal_vs_Acinar_core_paired.csv")

corecols <- InSituType::colorCellTypes(names = unique(psb$meta_data$tma_core_number), palette = "brewers")
coreidx <- psb$meta_data |> arrange(group, tma_core_number) |> rownames()
toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- (log2(DESeq2::counts(dds, normalized = T)+1)[toplot, coreidx] |> apply(MARGIN = 1, FUN = scale) |> t())
colnames(mat) <- psb$meta_data[coreidx,]$roi
ha <- HeatmapAnnotation(sub_type = psb$meta_data[coreidx,]$group,
                        tma_core = psb$meta_data[coreidx,]$tma_core_number,
                        col = list("sub_type"=c("Acinar"="orange", "Ductal"="yellow2"), 
                                   "tma_core"=corecols))
pdf(file = "10.4 - DESeq2_Ductal_vs_Acinar_core_paired_top_genes_heatmap.pdf", width = 8, height = 14)
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

# * Note that the batch corrected expression will require some extra work and foresight to get *

## 11) Intraductal_spread vs. Precursor ----------------------------------------
meta$sub_types_v5 <- case_when(meta$roi %in% c(48, 49) ~ "Precursor", 
                                        meta$sub_types_v2 %in% c("Acinar IDC-P", "Acinar IDC-P_crib", "AIP") ~ "Intraductal_spread", 
                                        T ~ meta$sub_types_v2)

small_meta <- meta |> filter(sub_types_v5 %in% c("Precursor", "Intraductal_spread"))
small_meta$group <- small_meta$sub_types_v5

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(small_meta)], colData = small_meta, design = ~group)
dds@colData$sizeFactor <- small_meta$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Precursor", ] |> colMeans()) - (mm[dds@colData$group == "Intraductal_spread", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 158
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.843003

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
  annotate(geom = "text", x = 0.6, y = -2, label = "Precursor", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Intraductal_spread", hjust = 1, fontface = "bold")
ggsave(filename = "11 - DESeq2_Precursor_vs_Intraductal_spread.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "11 - DESeq2_Precursor_vs_Intraductal_spread.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
idx <- sprintf("%03d", small_meta$roi)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = small_meta$group,
                        col = list("sub_type"=c("Precursor"="dodgerblue", "Intraductal_spread"="gold3")))
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

# * Note this comparison is INVALID! It is perfectly confounded by patient. *

## 12) Ductal vs. Ductal IDC-P -------------------------------------------------
idx <- meta$sub_types_v2 %in% c("Ductal", "Ductal IDC-P")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "Ductal", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "Ductal IDC-P", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 65
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.006423

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Ductal IDC-P", hjust = 1, fontface = "bold")
ggsave(filename = "12 - DESeq2_Ductal_vs_Ductal_IDC-P.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "12 - DESeq2_Ductal_vs_Ductal_IDC-P.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("Ductal"="yellow2", "Ductal IDC-P"="lightgreen")))
pdf(file = "12 - DESeq2_Ductal_vs_Ductal_IDC-P_top_genes_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
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
pdf(file = "12 - DESeq2_Ductal_vs_Ductal_IDC-P_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = meta[idx,]$sub_types_v2, top_annotation = ha,
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

## 13) Paired Ductal vs Ductal IDC-P -------------------------------------------
pairs <- table(meta$tma_core_number, meta$sub_types_v2) |> 
  as.data.frame.matrix() |> 
  select(Ductal, `Ductal IDC-P`)
pairs <- pairs[rowSums(pairs) > 1,]

paired_meta <- meta[meta$tma_core_number %in% rownames(pairs),] |> 
  filter(sub_types_v2 %in% c("Ductal", "Ductal IDC-P"))
paired_meta$group <- paired_meta$sub_types_v2

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,rownames(paired_meta)], colData = paired_meta, design = ~group+tma_core_number)
dds@colData$sizeFactor <- paired_meta$q_norm_qfactors

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~~group+tma_core_number, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "Ductal", ] |> colMeans()) - (mm[dds@colData$group == "Ductal IDC-P", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 0
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.045845

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
  annotate(geom = "text", x = -0.6, y = -2, label = "Ductal IDC-P", hjust = 1, fontface = "bold")
ggsave(filename = "13 - DESeq2_Ductal_vs_Ductal_IDC-P_core_paired.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "13 - DESeq2_Ductal_vs_Ductal_IDC-P_core_paired.csv")

## 14.1) All Acinar IDC-P vs all acinar carcinoma ------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  meta$sub_types_v2 %in% c("Acinar_non_crib_LG", "Acinar_non_crib_HG", "Acinar_crib") ~ "all_acinar_carcinoma",
  T ~ meta$sub_types_v2
)
idx <- meta$group %in% c("all_acinar_IDC-P", "all_acinar_carcinoma")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "all_acinar_carcinoma", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 87
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.031432

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
idx <- meta$group %in% c("all_acinar_IDC-P", "Acinar_non_crib_LG")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib_LG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 91
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9462617

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
idx <- meta$group %in% c("all_acinar_IDC-P", "Acinar_non_crib_HG")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_non_crib_HG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 205
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.574519

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
idx <- meta$group %in% c("all_acinar_IDC-P", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 13
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.143225

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

## 14.5) Acinar IDC-P_Crib vs. Acinar_Crib -------------------------------------
idx <- meta$sub_types %in% c("Acinar IDC-P_crib", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "Acinar IDC-P_crib", ] |> colMeans()) - (mm[dds@colData$sub_types == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 0
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9500585

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
  annotate(geom = "text", x = 0.6, y = -2, label = "Acinar IDC-P", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "14.5 - DESeq2_Acinar_IDC-P_crib_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "14.5 - DESeq2_Acinar_IDC-P_crib_vs_Acinar_crib.csv")

## 15.1) AIP vs. all acinar carcinoma ------------------------------------------
idx <- meta$sub_types_v2 %in% c("AIP", "Acinar_non_crib_LG", "Acinar_non_crib_HG", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types_v2+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$sub_types_v2 == "AIP", ] |> colMeans()) - 
  ((mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans())+(mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans())+(mm[dds@colData$sub_types_v2 == "Acinar_crib", ] |> colMeans()))/3
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 132
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.526369

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
  annotate(geom = "text", x = 0.6, y = -2, label = "AIP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "all_acinar_carcinoma", hjust = 1, fontface = "bold")
ggsave(filename = "15.1 - DESeq2_AIP_vs_all_acinar_carcinoma.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "15.1 - DESeq2_AIP_vs_all_acinar_carcinoma.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("AIP"="purple", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "15.1 - DESeq2_AIP_vs_all_acinar_carcinoma_top_genes_heatmap.pdf", width = 8, height = 14)
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
pdf(file = "15.1 - DESeq2_AIP_vs_all_acinar_carcinoma_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
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

## 15.2) AIP vs. all Acinar_non-crib_LG ----------------------------------------
idx <- meta$sub_types_v2 %in% c("AIP", "Acinar_non_crib_LG")
contr <- 
  (mm[dds@colData$sub_types_v2 == "AIP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_LG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 133
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.177462

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
  annotate(geom = "text", x = 0.6, y = -2, label = "AIP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_LG", hjust = 1, fontface = "bold")
ggsave(filename = "15.2 - DESeq2_AIP_vs_Acinar_non_crib_LG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "15.2 - DESeq2_AIP_vs_Acinar_non_crib_LG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("AIP"="purple", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "15.2 - DESeq2_AIP_vs_Acinar_non_crib_LG_top_genes_heatmap.pdf", width = 8, height = 14)
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
pdf(file = "15.2 - DESeq2_AIP_vs_Acinar_non_crib_LG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
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

## 15.3) AIP vs. all Acinar_non-crib_HG ----------------------------------------
idx <- meta$sub_types_v2 %in% c("AIP", "Acinar_non_crib_HG")
contr <- 
  (mm[dds@colData$sub_types_v2 == "AIP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "Acinar_non_crib_HG", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 430
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.768464

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
  annotate(geom = "text", x = 0.6, y = -2, label = "AIP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_non_crib_HG", hjust = 1, fontface = "bold")
ggsave(filename = "15.3 - DESeq2_AIP_vs_Acinar_non_crib_HG.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "15.3 - DESeq2_AIP_vs_Acinar_non_crib_HG.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("AIP"="purple", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "15.3 - DESeq2_AIP_vs_Acinar_non_crib_HG_top_genes_heatmap.pdf", width = 8, height = 14)
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
pdf(file = "15.3 - DESeq2_AIP_vs_Acinar_non_crib_HG_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
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

## 15.4) AIP vs. all Acinar_crib -----------------------------------------------
idx <- meta$sub_types_v2 %in% c("AIP", "Acinar_crib")
contr <- 
  (mm[dds@colData$sub_types_v2 == "AIP", ] |> colMeans()) - (mm[dds@colData$sub_types_v2 == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 54
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 1.399424

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
  annotate(geom = "text", x = 0.6, y = -2, label = "AIP", hjust = 0, fontface = "bold") +
  annotate(geom = "segment", x = 0, xend = -0.5, y = -2, yend = -2, arrow = arrow(type = "open", length = unit(4, "pt"))) +
  annotate(geom = "text", x = -0.6, y = -2, label = "Acinar_crib", hjust = 1, fontface = "bold")
ggsave(filename = "15.4 - DESeq2_AIP_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
write.csv(x = res, file = "15.4 - DESeq2_AIP_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types_v2,
                        col = list("sub_type"=c("AIP"="purple", "Acinar_non_crib_LG"="pink", "Acinar_non_crib_HG"="red3", "Acinar_crib"="green4")))
pdf(file = "15.4 - DESeq2_AIP_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
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
pdf(file = "15.4 - DESeq2_AIP_vs_Acinar_crib_top_genes_heatmap_batch_corrected.pdf", width = 8, height = 14)
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

## 16) All Acinar IDC-P vs AIP -------------------------------------------------
meta$group <- case_when(
  meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "Acinar IDC-P") ~ "all_acinar_IDC-P",
  T ~ meta$sub_types_v2
)
idx <- meta$group %in% c("all_acinar_IDC-P", "AIP")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~group+patient_deid)

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~group+patient_deid, data = dds@colData)
contr <- 
  (mm[dds@colData$group == "all_acinar_IDC-P", ] |> colMeans()) - (mm[dds@colData$group == "AIP", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 52
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.7422728

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

## Session ---------------------------------------------------------------------
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
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
#   [1] presto_1.0.0          data.table_1.16.4     Rcpp_1.0.14           lme4_1.1-37           Matrix_1.7-3          patchwork_1.3.0      
# [7] magrittr_2.0.3        ComplexHeatmap_2.22.0 ggplot2_3.5.2         dplyr_1.1.4          
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0              shape_1.4.6.1              
# [5] umap_0.2.10.0               spatstat.utils_3.1-3        ggbeeswarm_0.7.2            farver_2.1.2               
# [9] nloptr_2.2.1                GlobalOptions_0.1.2         zlibbioc_1.52.0             vctrs_0.6.5                
# [13] minqa_1.2.8                 askpass_1.2.1               htmltools_0.5.8.1           S4Arrays_1.6.0             
# [17] cellranger_1.1.0            SparseArray_1.6.2           parallelly_1.43.0           NanoStringNCTools_1.14.0   
# [21] htmlwidgets_1.6.4           plyr_1.8.9                  uuid_1.2-1                  lifecycle_1.0.4            
# [25] iterators_1.0.14            pkgconfig_2.0.3             R6_2.6.1                    fastmap_1.2.0              
# [29] future_1.40.0               GenomeInfoDbData_1.2.13     rbibutils_2.3               MatrixGenerics_1.18.1      
# [33] clue_0.3-66                 numDeriv_2016.8-1.1         digest_0.6.37               ggnewscale_0.5.1           
# [37] GGally_2.2.1                colorspace_2.1-1            S4Vectors_0.44.0            DESeq2_1.46.0              
# [41] RSpectra_0.16-2             irlba_2.3.5.1               InSituType_2.0              GenomicRanges_1.58.0       
# [45] SnowballC_0.7.1             labeling_0.4.3              progressr_0.15.1            httr_1.4.7                 
# [49] polyclip_1.10-7             abind_1.4-8                 compiler_4.4.2              withr_3.0.2                
# [53] doParallel_1.0.17           BiocParallel_1.40.2         viridis_0.6.5               ggstats_0.9.0              
# [57] MASS_7.3-61                 openssl_2.3.2               DelayedArray_0.32.0         rjson_0.2.23               
# [61] tools_4.4.2                 vipor_0.4.7                 beeswarm_0.4.0              future.apply_1.11.3        
# [65] glue_1.8.0                  nlme_3.1-166                cluster_2.1.6               reshape2_1.4.4             
# [69] generics_0.1.3              gtable_0.3.6                spatstat.data_3.1-6         tidyr_1.3.1                
# [73] sp_2.2-0                    XVector_0.46.0              BiocGenerics_0.52.0         spatstat.geom_3.3-6        
# [77] ggrepel_0.9.6               foreach_1.5.2               pillar_1.10.1               stringr_1.5.1              
# [81] spam_2.11-1                 limma_3.62.2                circlize_0.4.16             splines_4.4.2              
# [85] lattice_0.22-6              renv_1.1.1                  deldir_2.0-4                tidyselect_1.2.1           
# [89] SingleCellExperiment_1.28.1 locfit_1.5-9.12             Biostrings_2.74.1           reformulas_0.4.0           
# [93] gridExtra_2.3               IRanges_2.40.1              SummarizedExperiment_1.36.0 scattermore_1.2            
# [97] stats4_4.4.2                Biobase_2.66.0              statmod_1.5.0               matrixStats_1.5.0          
# [101] pheatmap_1.0.12             stringi_1.8.4               UCSC.utils_1.2.0            boot_1.3-31                
# [105] codetools_0.2-20            lsa_0.73.3                  tibble_3.2.1                BiocManager_1.30.25        
# [109] cli_3.6.4                   uwot_0.2.3                  reticulate_1.42.0           systemfonts_1.2.3          
# [113] Rdpack_2.6.4                munsell_0.5.1               GenomeInfoDb_1.42.3         readxl_1.4.5               
# [117] globals_0.17.0              EnvStats_3.0.0              png_0.1-8                   spatstat.univar_3.1-2      
# [121] parallel_4.4.2              GeomxTools_3.10.0           dotCall64_1.2               mclust_6.1.1               
# [125] listenv_0.9.1               viridisLite_0.4.2           ggthemes_5.1.0              lmerTest_3.1-3             
# [129] ggiraph_0.8.13              scales_1.3.0                SeuratObject_5.0.2          purrr_1.0.4                
# [133] crayon_1.5.3                GetoptLong_1.0.5            rlang_1.1.5                

