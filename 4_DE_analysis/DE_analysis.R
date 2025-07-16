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
meta <- read.csv(file = "meta.csv", row.names = 1)

# Remove NegProbe-WTX before doing any of this analysis
norm <- norm[rownames(norm) != "NegProbe-WTX",]
cts <- cts[rownames(cts) != "NegProbe-WTX",]

## From Amaya ## ---------------------------------------------------------------
# For DE analysis, Lets start with these -
# Prostate cancer -
# 1. Ductal vs Acinar non-crib
# 2. Ductal vs Acinar crib
# 3. Ductal IDC-P vs Acinar IDC-P
# 4. Ductal vs Ductal IDC-P
# 5. Acinar non-crib vs crib
# ------------------------------------------------------------------------------

# We need the patient identifier for each AOI.
tma_ls <- list()
for (tma in c("TMA1", "TMA2")) {
  tma_ls[[tma]] <- openxlsx::read.xlsx(xlsxFile = "TMA_maps.xlsx", 
                                      sheet = tma, startRow = 15, 
                                      cols = 9:10)
}
samplevelmeta <- dplyr::bind_rows(tma_ls, .id = "tma")
samplevelmeta$patient_deid <- paste("pt", as.factor(samplevelmeta$patient) |> as.numeric(), sep = "_")
# write.csv(x = samplevelmeta, file = "patient_map.csv")
# write.csv(x = samplevelmeta[,c("tma", "core", "patient_deid")], file = "patient_map_de-identified.csv")

# Adding patient IDs to the data
samplevelmeta$tmacore <- paste(samplevelmeta$tma, samplevelmeta$core, sep = "")
meta$patient_deid <- plyr::mapvalues(x = meta$tma_core_number, from = samplevelmeta$tmacore, to = samplevelmeta$patient_deid)
# write.csv(x = meta, file = "meta.csv")

# Let's look at the study design
table(meta$sub_types, meta$patient_deid)
table(meta$sub_types, meta$cancer_type)


# Each patient has 1-4 observations. Most of them have 2.
table(meta$sub_types, meta$patient_deid) |> colSums() |> table()
# 1  2  3  4 
# 6 12  6 10 

# Of those that have 3 or 4, most of the patients have each sampled AOI from different categories. 
# This is not ideal experimental design.
# OPTIONS:
# 1) Mixed model --> probably not good, since many patients are represented by 1 observation per grouping
# 2) Pseudo-bulk --> this would be ideal for when either all AOIs from each patient belong to the same group or 
#                    for when each patient has observations in the groups being compared --> probably not good here
# 3) Fixed effects model with patient as a predictor. Probably okay, but it will not really account for all the bias.

# Originally, I thought that it made sense to pseudo-bulk what we could, then do a fixed effects model with patient as a predictor. 
# This became complicated quickly, so I just went with option 3: 
# DESeq2 with this formula: ~subtype+patient

# Modeling and testing:
## Ductal vs Acinar non-crib ---------------------------------------------------
meta$sub_types <- gsub(pattern = "\\-", replacement = "_", x = meta$sub_types)
idx <- meta$sub_types %in% c("Ductal", "Acinar_non_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types+patient_deid)

# We could do this, but I want to use the exact same normalization factors as we used to produce our normalized data.
# nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
# sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
# sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
# dds@colData$sizeFactor <- sfs

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types+patient_deid, data = dds@colData)
contr <- (mm[dds@colData$sub_types == "Ductal", ] |> colMeans()) - 
  (mm[dds@colData$sub_types == "Acinar_non_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame() # Note: set independentFiltering = F to prevent NA values if needed
res$padj |> is.na() |> mean()
# [1] 0.2132853
sum(res$padj < 0.05, na.rm = T)
# [1] 515
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.8848523

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
# ggsave(filename = "DESeq2_Ductal_vs_Acinar_non_crib.pdf", device = "pdf", width = 8, height = 8)
# write.csv(x = res, file = "DESeq2_Ductal_vs_Acinar_non_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- (norm[toplot, idx] |> apply(MARGIN = 1, FUN = scale) |> t())
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types, col = list("sub_type"=c("Ductal"="gold", "Acinar_non_crib"="red3")))
# pdf(file = "DESeq2_Ductal_vs_Acinar_non_crib_top_genes_heatmap.pdf", width = 8, height = 14)
ComplexHeatmap::Heatmap(matrix = mat, column_split = meta[idx,]$sub_types, top_annotation = ha,
                        row_names_gp = gpar(fontsize = 7), show_column_names = F,
                        name = "row-scaled\nlog2(q3norm+1)",
                        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.02), quantile(mat, 0.98), length.out=51),
                                                   colors = viridis::viridis(51)), 
                        width = ncol(mat)*unit(3, "mm"),
                        height = nrow(mat)*unit(2, "mm")
)
# dev.off()

## Ductal vs Acinar crib ---------------------------------------------------
idx <- meta$sub_types %in% c("Ductal", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types+patient_deid)

# nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
# sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
# sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
# dds@colData$sizeFactor <- sfs

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types+patient_deid, data = dds@colData)
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
# ggsave(filename = "DESeq2_Ductal_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
# write.csv(x = res, file = "DESeq2_Ductal_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- (norm[toplot, idx] |> apply(MARGIN = 1, FUN = scale) |> t())
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types, col = list("sub_type"=c("Ductal"="gold", "Acinar_crib"="green4")))
# pdf(file = "DESeq2_Ductal_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
ComplexHeatmap::Heatmap(matrix = mat, column_split = meta[idx,]$sub_types, top_annotation = ha,
                        row_names_gp = gpar(fontsize = 7), show_column_names = F,
                        name = "row-scaled\nlog2(q3norm+1)",
                        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.02), quantile(mat, 0.98), length.out=51),
                                                   colors = viridis::viridis(51)), 
                        width = ncol(mat)*unit(3, "mm"),
                        height = nrow(mat)*unit(2, "mm")
)
# dev.off()

## Acinar non-crib vs Acinar crib ---------------------------------------------------
idx <- meta$sub_types %in% c("Acinar_non_crib", "Acinar_crib")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts[,idx], colData = meta[idx,], design = ~sub_types+patient_deid)

# nfs <- dds@assays@data$counts |> edgeR::DGEList() |> edgeR::calcNormFactors(method = "upperquartile")
# sfs <- (nfs$samples$lib.size*nfs$samples$norm.factors)
# sfs <- sfs / exp(mean(log(sfs))) # Scaling by geometric mean, which DESeq2 expects
# dds@colData$sizeFactor <- sfs

dds@colData$sizeFactor <- meta$q_norm_qfactors[idx]

dds <- DESeq2::DESeq(object = dds)
DESeq2::plotDispEsts(dds)
mm <- model.matrix(~sub_types+patient_deid, data = dds@colData)
contr <- (mm[dds@colData$sub_types == "Acinar_non_crib", ] |> colMeans()) - 
  (mm[dds@colData$sub_types == "Acinar_crib", ] |> colMeans())
res <- DESeq2::results(object = dds, contrast = contr) |> as.data.frame()
sum(res$padj < 0.05, na.rm = T)
# [1] 151
sum(res$log2FoldChange < 0)/sum(res$log2FoldChange > 0)
# [1] 0.9941482

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
# ggsave(filename = "DESeq2_Acinar_non_crib_vs_Acinar_crib.pdf", device = "pdf", width = 8, height = 8)
# write.csv(x = res, file = "DESeq2_Acinar_non_crib_vs_Acinar_crib.csv")

toplot <- na.omit(res) |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) |> pull(target)
mat <- (norm[toplot, idx] |> apply(MARGIN = 1, FUN = scale) |> t())
ha <- HeatmapAnnotation(sub_type = meta[idx,]$sub_types, col = list("sub_type"=c("Acinar_non_crib"="red3", "Acinar_crib"="green4")))
# pdf(file = "DESeq2_Acinar_non_crib_vs_Acinar_crib_top_genes_heatmap.pdf", width = 8, height = 14)
ComplexHeatmap::Heatmap(matrix = mat, column_split = meta[idx,]$sub_types, top_annotation = ha,
                        row_names_gp = gpar(fontsize = 7), show_column_names = F,
                        name = "row-scaled\nlog2(q3norm+1)",
                        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.02), quantile(mat, 0.98), length.out=51),
                                                   colors = viridis::viridis(51)), 
                        width = ncol(mat)*unit(3, "mm"),
                        height = nrow(mat)*unit(2, "mm")
)
# dev.off()

# Session
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
#   [1] tidyselect_1.2.1            viridisLite_0.4.2           farver_2.1.2                viridis_0.6.5              
# [5] digest_0.6.37               lifecycle_1.0.4             cluster_2.1.6               statmod_1.5.0              
# [9] compiler_4.4.2              rlang_1.1.5                 tools_4.4.2                 labeling_0.4.3             
# [13] S4Arrays_1.6.0              DelayedArray_0.32.0         plyr_1.8.9                  RColorBrewer_1.1-3         
# [17] abind_1.4-8                 BiocParallel_1.40.2         withr_3.0.2                 purrr_1.0.4                
# [21] BiocGenerics_0.52.0         stats4_4.4.2                colorspace_2.1-1            edgeR_4.4.2                
# [25] scales_1.3.0                iterators_1.0.14            MASS_7.3-61                 SummarizedExperiment_1.36.0
# [29] cli_3.6.4                   crayon_1.5.3                reformulas_0.4.0            generics_0.1.3             
# [33] httr_1.4.7                  rjson_0.2.23                minqa_1.2.8                 stringr_1.5.1              
# [37] zlibbioc_1.52.0             splines_4.4.2               ggthemes_5.1.0              parallel_4.4.2             
# [41] BiocManager_1.30.25         XVector_0.46.0              matrixStats_1.5.0           vctrs_0.6.5                
# [45] boot_1.3-31                 jsonlite_2.0.0              IRanges_2.40.1              GetoptLong_1.0.5           
# [49] S4Vectors_0.44.0            ggrepel_0.9.6               clue_0.3-66                 scattermore_1.2            
# [53] locfit_1.5-9.12             foreach_1.5.2               limma_3.62.2                glue_1.8.0                 
# [57] nloptr_2.2.1                codetools_0.2-20            stringi_1.8.4               shape_1.4.6.1              
# [61] gtable_0.3.6                GenomeInfoDb_1.42.3         GenomicRanges_1.58.0        UCSC.utils_1.2.0           
# [65] munsell_0.5.1               tibble_3.2.1                pillar_1.10.1               GenomeInfoDbData_1.2.13    
# [69] circlize_0.4.16             R6_2.6.1                    Rdpack_2.6.4                doParallel_1.0.17          
# [73] lattice_0.22-6              Biobase_2.66.0              rbibutils_2.3               png_0.1-8                  
# [77] openxlsx_4.2.8              renv_1.1.1                  zip_2.3.2                   gridExtra_2.3              
# [81] SparseArray_1.6.2           nlme_3.1-166                DESeq2_1.46.0               MatrixGenerics_1.18.1      
# [85] pkgconfig_2.0.3             GlobalOptions_0.1.2        

