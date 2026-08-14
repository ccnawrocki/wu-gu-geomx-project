rm(list = ls())

# Bladder Data
cts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("meta_amended.csv", row.names = 1)
classres <- openxlsx::read.xlsx(xlsxFile = "5_tissue_classification/consensusMIBC_classification_results.xlsx", sheet = 1, rowNames = T)

##### !!! #####
# Need to remove the normal AOI!
bmeta <- meta[meta$cancer_type == "bladder",]
bmeta$MolecularClass <- plyr::mapvalues(x = bmeta$sample_id, from = rownames(classres), to = classres$consensusClass)
bmeta <- bmeta[bmeta$sub_types_v2 != "Normal urothelium",]

# Top 10 markers for each molecular class
library(limma)
library(edgeR)
bcts <- cts[rownames(cts) != "NegProbe-WTX" ,rownames(bmeta)]
y <- DGEList(counts = bcts) |> calcNormFactors(method = "upperquartile")
mm <- model.matrix(~MolecularClass, data = bmeta)
v <- voomLmFit(counts = y, 
               design = mm, 
               sample.weights = T,
               plot = T, 
               save.plot = T, 
               keep.EList = T
)

listofparams <- list()
for (MC in unique(bmeta$MolecularClass)) {
  listofparams[[MC]] <- mm[bmeta$MolecularClass == MC,] |> colMeans()
}
parameterizations <- listofparams |> dplyr::bind_rows() |> as.matrix()
rownames(parameterizations) <- names(listofparams)

outs <- list()
for (MC in unique(bmeta$MolecularClass)) {
  MCoi <- parameterizations[MC,]
  MCnoi <- (parameterizations[rownames(parameterizations) != MC,] |> colSums())/5
  contr <- (MCoi-MCnoi)
  out <- contrasts.fit(fit = v, contrasts = contr) |> eBayes() |> topTable(n = Inf)
  outs[[MC]] <- out
}

outs <- lapply(X = outs, FUN = function(df) {df$target <- rownames(df); return(df)})
markers <- dplyr::bind_rows(outs, .id = "MolecularClass") |> magrittr::set_rownames(NULL)

## - Need to have p<0.05
## - Need to have logFC>1

library(dplyr)
topmarks <- filter(markers, logFC > 0) |> 
  mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH")) |> 
  filter(adj.P.Val < 0.01 & logFC > 1) |> 
  group_by(MolecularClass) |> 
  top_n(n = -10, wt = adj.P.Val) |> 
  arrange(MolecularClass, desc(logFC))
#write.csv(x = topmarks, file = "10_revisions/bladder/MolecularClass_Top10Markers_updated_limma.csv", row.names = F)

# Visualization
mindray_discrete <- c("#10BBB9", "#233987", "#030306", "#C7360B", "#F6BF15")
mindray_palette <- colorRampPalette(mindray_discrete)(256)
image(matrix(1:256, ncol=1), col = mindray_palette, axes = FALSE, main = "Mindray Doppler Palette")

library(ComplexHeatmap)

mat <- v$EList$E[topmarks$target, bmeta |> arrange(MolecularClass, patient_deid) |> pull(sample_id)] |> apply(MARGIN = 1, FUN = scale) |> t()
bmeta$sub_types_v2 <- ifelse(test = bmeta$sub_types_v2 == "small cell", yes = "SCNEC", no = bmeta$sub_types_v2)

pdf(file = "10_revisions/bladder/MolecularClass_Top10Markers_with_patient_updated_limma.pdf", width = 8, height = 10)
Heatmap(mat, 
        top_annotation = HeatmapAnnotation(Class = bmeta |> arrange(MolecularClass, patient_deid) |> pull(MolecularClass), 
                                           Pathology = bmeta |> arrange(MolecularClass, patient_deid) |> pull(sub_types_v2),
                                           Patient = bmeta |> arrange(MolecularClass, patient_deid) |> pull(patient_deid),
                                           col = list("Class" = 
                                                        c("Ba/Sq" = "pink", 
                                                          "LumNS" = "orange", 
                                                          "LumP" = "gold", 
                                                          "LumU" = "limegreen", 
                                                          "NE-like" = "grey", 
                                                          "Stroma-rich" = "dodgerblue"), 
                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(bmeta$patient_deid)] |> 
                                                        setNames(nm = unique(bmeta$patient_deid)), 
                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue")
                                           ),
                                           show_annotation_name = F),
        cluster_columns = F, cluster_rows = F,
        col = circlize::colorRamp2(breaks = seq(-4, 4, length.out = 101), 
                                   colors = colorRampPalette(mindray_discrete)(101)), 
        name = "Scaled\nExpression")
dev.off()

# What if I use DESeq2 to keep the methods consistent? 
dds <- DESeq2::DESeqDataSetFromMatrix(countData = bcts, 
                                      colData = bmeta, 
                                      design = mm)
DESeq2::sizeFactors(dds) <- dds$q_norm_qfactors
dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)

outs <- list()
for (MC in unique(bmeta$MolecularClass)) {
  MCoi <- parameterizations[MC,]
  MCnoi <- (parameterizations[rownames(parameterizations) != MC,] |> colSums())/5
  contr <- (MCoi-MCnoi)
  out <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
  outs[[MC]] <- out
}
outs <- lapply(X = outs, FUN = function(df) {df$target <- rownames(df); return(df)})
markers <- dplyr::bind_rows(outs, .id = "MolecularClass") |> magrittr::set_rownames(NULL)

topmarks <- filter(markers, log2FoldChange > 0) |> 
  mutate(adj.P.Val = p.adjust(p = pvalue, method = "BH")) |> 
  filter(adj.P.Val < 0.01 & log2FoldChange > 1) |> 
  group_by(MolecularClass) |> 
  top_n(n = -10, wt = adj.P.Val) |> 
  arrange(MolecularClass, desc(log2FoldChange))
#write.csv(x = topmarks, file = "10_revisions/bladder/MolecularClass_Top10Markers_updated_deseq2.csv", row.names = F)

norm <- read.csv(file = "log2plus1_q3norm.csv", row.names = 1)
mat <- norm[topmarks$target, bmeta |> arrange(MolecularClass, patient_deid) |> pull(sample_id)] |> apply(MARGIN = 1, FUN = scale) |> t()
bmeta$sub_types_v2 <- ifelse(test = bmeta$sub_types_v2 == "small cell", yes = "SCNEC", no = bmeta$sub_types_v2)

pdf(file = "10_revisions/bladder/MolecularClass_Top10Markers_with_patient_updated_deseq2.pdf", width = 8, height = 10)
Heatmap(mat, 
        top_annotation = HeatmapAnnotation(Class = bmeta |> arrange(MolecularClass, patient_deid) |> pull(MolecularClass), 
                                           Pathology = bmeta |> arrange(MolecularClass, patient_deid) |> pull(sub_types_v2),
                                           Patient = bmeta |> arrange(MolecularClass, patient_deid) |> pull(patient_deid),
                                           col = list("Class" = 
                                                        c("Ba/Sq" = "pink", 
                                                          "LumNS" = "orange", 
                                                          "LumP" = "gold", 
                                                          "LumU" = "limegreen", 
                                                          "NE-like" = "grey", 
                                                          "Stroma-rich" = "dodgerblue"), 
                                                      "Patient" = ggprism::ggprism_data$colour_palettes$colors[1:dplyr::n_distinct(bmeta$patient_deid)] |> 
                                                        setNames(nm = unique(bmeta$patient_deid)), 
                                                      "Pathology" = c("CIUC" = "pink3", "UCIS" = "cyan", "MP" = "yellow2", "PC" = "grey", "PUC" = "limegreen", "SM" = "black", "SCNEC" = "red", "Normal" = "blue")
                                           ),
                                           show_annotation_name = T),
        cluster_columns = F, cluster_rows = F,
        col = circlize::colorRamp2(breaks = seq(-4, 4, length.out = 101), 
                                   colors = colorRampPalette(mindray_discrete)(101)), 
        name = "Scaled\nExpression")
dev.off()

