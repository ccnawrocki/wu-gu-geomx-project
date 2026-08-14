rm(list = ls())
.rs.restartR()

cts <- read.csv(file = "counts.csv", row.names = 1)
meta <- read.csv(file = "meta_amended.csv", row.names = 1)

classification <- openxlsx::read.xlsx(xlsxFile = "5_tissue_classification/consensusMIBC_classification_results.xlsx", rowNames = T)
meta$class <- plyr::mapvalues(x = rownames(meta), from = rownames(classification), to = classification$consensusClass)

bmeta <- meta[meta$cancer_type == "bladder" & meta$sub_types_v2 != "Normal urothelium",]
bcts <- cts[rownames(cts) != "NegProbe-WTX", rownames(bmeta)]


## START WITH OVERALL CLUSTERS ## ----------------------------------------------
mm <- model.matrix(~1, data = bmeta)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = bcts, 
                                      colData = bmeta, 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

hvgs <- (dds@rowRanges@elementMetadata@listData |> as.data.frame() |> dplyr::arrange(desc(dispGeneEst)) |> rownames())[1:2500]

norm <- read.csv(file = "log2plus1_q3norm.csv", row.names = 1)
mat <- norm[hvgs, rownames(bmeta)] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", bmeta$roi)

clusts <- hclust(d = dist(t(mat)))
clusts <- cutree(tree = clusts, k = 3)

mindray_discrete <- c("#10BBB9", "#233987", "#030306", "#C7360B", "#F6BF15")
mindray_palette <- colorRampPalette(mindray_discrete)(256)

bmeta$sub_types_v2 <- ifelse(test = bmeta$sub_types_v2 == "small cell", yes = "SCNEC", no = bmeta$sub_types_v2)

library(ComplexHeatmap)

pdf(file = "10_revisions/bladder/all_rois_unsupervised_clusters_from_2500_hvgs.pdf", width = 6, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts,
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = F, column_names_gp = gpar(fontsize = 6),
                        top_annotation = HeatmapAnnotation(Class = bmeta$class, 
                                                           Pathology = bmeta$sub_types_v2,
                                                           Patient = bmeta$patient_deid, 
                                                           Segment = bmeta$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(bmeta$patient_deid)] |> setNames(nm = unique(bmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
                        )
dev.off()

bmeta$cluster <- clusts |> as.character()
mm <- model.matrix(~cluster, data = bmeta[bmeta$cluster %in% 1:2,])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = bcts[,bmeta$cluster %in% 1:2], 
                                      colData = bmeta[bmeta$cluster %in% 1:2,], 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

markers <-  DESeq2::results(object = dds, name = "cluster2") |> as.data.frame()
markers$target <- rownames(markers)
topmarks <- filter(markers, abs(log2FoldChange) > 1) |> 
  mutate(enrichment = ifelse(test = log2FoldChange > 0, yes = "cluster2", no = "cluster1")) |> 
  filter(padj < 0.01) |> 
  group_by(enrichment) |> 
  top_n(n = -25, wt = padj) |> 
  arrange(enrichment, desc(log2FoldChange))

mat <- norm[topmarks$target, rownames(bmeta)][,bmeta$cluster %in% 1:2] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", bmeta[bmeta$cluster %in% 1:2,]$roi)

pdf(file = "10_revisions/bladder/all_rois_clusters_1_and_2_top25_markers.pdf", width = 6, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts[bmeta$cluster %in% 1:2],
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = T, column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = bmeta[bmeta$cluster %in% 1:2,]$class, 
                                                           Pathology = bmeta[bmeta$cluster %in% 1:2,]$sub_types_v2,
                                                           Patient = bmeta[bmeta$cluster %in% 1:2,]$patient_deid, 
                                                           Segment = bmeta[bmeta$cluster %in% 1:2,]$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(bmeta$patient_deid)] |> setNames(nm = unique(bmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()


## NEXT, JUST STROMAL DOMINANT ## ----------------------------------------------
fbmeta <- bmeta[bmeta$cluster == "1",]
fbcts <- bcts[,rownames(fbmeta)]
mm <- model.matrix(~1, data = fbmeta)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = fbcts, 
                                      colData = fbmeta, 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

hvgs <- (dds@rowRanges@elementMetadata@listData |> as.data.frame() |> dplyr::arrange(desc(dispGeneEst)) |> rownames())[1:2500]

norm <- read.csv(file = "log2plus1_q3norm.csv", row.names = 1)
mat <- norm[hvgs, rownames(fbmeta)] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta$roi)

clusts <- hclust(d = dist(t(mat)))
clusts <- cutree(tree = clusts, k = 3)

pdf(file = "10_revisions/bladder/c1_stromal_dominant_unsupervised_subclusters_from_2500_hvgs.pdf", width = 4, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts,
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = F, column_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta$class, 
                                                           Pathology = fbmeta$sub_types_v2,
                                                           Patient = fbmeta$patient_deid, 
                                                           Segment = fbmeta$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()

fbmeta$subcluster <- clusts |> as.character()
mm <- model.matrix(~subcluster, data = fbmeta[fbmeta$subcluster %in% 1:2,])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = fbcts[,fbmeta$subcluster %in% 1:2], 
                                      colData = fbmeta[fbmeta$subcluster %in% 1:2,], 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

markers <-  DESeq2::results(object = dds, name = "subcluster2") |> as.data.frame()
markers$target <- rownames(markers)
topmarks <- filter(markers, abs(log2FoldChange) > 1) |> 
  mutate(enrichment = ifelse(test = log2FoldChange > 0, yes = "subcluster2", no = "subcluster1")) |> 
  filter(padj < 0.01) |> 
  group_by(enrichment) |> 
  top_n(n = -25, wt = padj) |> 
  arrange(enrichment, desc(log2FoldChange))

mat <- norm[topmarks$target, rownames(fbmeta)][,fbmeta$subcluster %in% 1:2] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta[fbmeta$subcluster %in% 1:2,]$roi)

pdf(file = "10_revisions/bladder/c1_stromal_dominant_clusters_1_and_2_markers.pdf", width = 8, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts[fbmeta$subcluster %in% 1:2],
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = T, column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta[fbmeta$subcluster %in% 1:2,]$class, 
                                                           Pathology = fbmeta[fbmeta$subcluster %in% 1:2,]$sub_types_v2,
                                                           Patient = fbmeta[fbmeta$subcluster %in% 1:2,]$patient_deid, 
                                                           Segment = fbmeta[fbmeta$subcluster %in% 1:2,]$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()

# Selected from MDA paper:
selected <- c("CD44", "KRT5", "KRT6A", "KRT14", "CDH3", "CD24", "FOXA1", "GATA3", "ERBB2", "ERBB3", "XBP1", "KRT18", "KRT19", "KRT20")

mat <- norm[selected, rownames(fbmeta)][,fbmeta$subcluster %in% 1:2] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta[fbmeta$subcluster %in% 1:2,]$roi)

pdf(file = "10_revisions/bladder/c1_stromal_dominant_clusters_1_and_2_MDA_selected_genes.pdf", width = 6, height = 6)
ComplexHeatmap::Heatmap(mat, column_split = clusts[fbmeta$subcluster %in% 1:2],
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = T, column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta[fbmeta$subcluster %in% 1:2,]$class, 
                                                           Pathology = fbmeta[fbmeta$subcluster %in% 1:2,]$sub_types_v2,
                                                           Patient = fbmeta[fbmeta$subcluster %in% 1:2,]$patient_deid, 
                                                           Segment = fbmeta[fbmeta$subcluster %in% 1:2,]$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()


## NEXT, PURE EPITHELIUM ## ----------------------------------------------------
fbmeta <- bmeta[bmeta$cluster == "2",]
fbcts <- bcts[,rownames(fbmeta)]
mm <- model.matrix(~1, data = fbmeta)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = fbcts, 
                                      colData = fbmeta, 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

hvgs <- (dds@rowRanges@elementMetadata@listData |> as.data.frame() |> dplyr::arrange(desc(dispGeneEst)) |> rownames())[1:2500]

norm <- read.csv(file = "log2plus1_q3norm.csv", row.names = 1)
mat <- norm[hvgs, rownames(fbmeta)] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta$roi)

clusts <- hclust(d = dist(t(mat)))
clusts <- cutree(tree = clusts, k = 2)

pdf(file = "10_revisions/bladder/c2_epithelium_pure_unsupervised_subclusters_from_2500_hvgs.pdf", width = 4, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts,
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = F, column_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta$class, 
                                                           Pathology = fbmeta$sub_types_v2,
                                                           Patient = fbmeta$patient_deid, 
                                                           Segment = fbmeta$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()

fbmeta$subcluster <- clusts |> as.character()
mm <- model.matrix(~subcluster, data = fbmeta[fbmeta$subcluster %in% 1:2,])
dds <- DESeq2::DESeqDataSetFromMatrix(countData = fbcts[,fbmeta$subcluster %in% 1:2], 
                                      colData = fbmeta[fbmeta$subcluster %in% 1:2,], 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

markers <-  DESeq2::results(object = dds, name = "subcluster2") |> as.data.frame()
markers$target <- rownames(markers)
topmarks <- filter(markers, abs(log2FoldChange) > 1) |> 
  mutate(enrichment = ifelse(test = log2FoldChange > 0, yes = "subcluster2", no = "subcluster1")) |> 
  filter(padj < 0.01) |> 
  group_by(enrichment) |> 
  top_n(n = -25, wt = padj) |> 
  arrange(enrichment, desc(log2FoldChange))

mat <- norm[topmarks$target, rownames(fbmeta)][,fbmeta$subcluster %in% 1:2] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta[fbmeta$subcluster %in% 1:2,]$roi)

pdf(file = "10_revisions/bladder/c2_epithelium_pure_clusters_1_and_2_markers.pdf", width = 8, height = 8)
ComplexHeatmap::Heatmap(mat, column_split = clusts[fbmeta$subcluster %in% 1:2],
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = T, column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta[fbmeta$subcluster %in% 1:2,]$class, 
                                                           Pathology = fbmeta[fbmeta$subcluster %in% 1:2,]$sub_types_v2,
                                                           Patient = fbmeta[fbmeta$subcluster %in% 1:2,]$patient_deid, 
                                                           Segment = fbmeta[fbmeta$subcluster %in% 1:2,]$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()

mat <- norm[selected, rownames(fbmeta)][,fbmeta$subcluster %in% 1:2] |> apply(MARGIN = 1, FUN = scale) |> t()
colnames(mat) <- sprintf(fmt = "%03s", fbmeta[fbmeta$subcluster %in% 1:2,]$roi)

pdf(file = "10_revisions/bladder/c2_epithelium_pure_clusters_1_and_2_MDA_selected_genes.pdf", width = 6, height = 6)
ComplexHeatmap::Heatmap(mat, column_split = clusts[fbmeta$subcluster %in% 1:2],
                        col = circlize::colorRamp2(breaks = seq(-3, 3, length.out = 71), colorRampPalette(mindray_discrete)(71)), 
                        name = "Scaled\nExpression", 
                        show_column_dend = T, show_row_dend = T, show_row_names = T, column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 8),
                        top_annotation = HeatmapAnnotation(Class = fbmeta[fbmeta$subcluster %in% 1:2,]$class, 
                                                           Pathology = fbmeta[fbmeta$subcluster %in% 1:2,]$sub_types_v2,
                                                           Patient = fbmeta[fbmeta$subcluster %in% 1:2,]$patient_deid, 
                                                           Segment = fbmeta[fbmeta$subcluster %in% 1:2,]$segment, 
                                                           col = list("Class" = c("Ba/Sq" = "pink", "LumNS" = "orange", "LumP" = "gold", "LumU" = "limegreen", "NE-like" = "grey", "Stroma-rich" = "dodgerblue"),
                                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue"), 
                                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(fbmeta$patient_deid)] |> setNames(nm = unique(fbmeta$patient_deid)), 
                                                                      "Segment" = c("Epithelium" = "grey", "Full ROI" = "black")))
)
dev.off()

