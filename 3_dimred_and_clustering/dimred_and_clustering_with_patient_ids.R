rm(list = ls())

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/wu-gu-prostate-bladder/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(ggbiplot)
library(magrittr)
library(patchwork)

# Data 
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta.csv", row.names = 1)

# Just want to check that I found size factors correctly
qs <- apply(X = cts, MARGIN = 2, FUN = quantile, 0.75)
nfs <- qs / EnvStats::geoMean(qs)
names(nfs) <- NULL
all.equal(nfs, meta$q_norm_qfactors)
# [1] TRUE

# Looks good. 
# Note: can also find geomean by doing exp(mean(log(qs)))
# Note: GeomxTools adds a threshold to the geomean under certain circumstances (not here though)

# Remove NegProbe-WTX before doing any of this analysis
norm <- norm[rownames(norm) != "NegProbe-WTX",]
cts <- cts[rownames(cts) != "NegProbe-WTX",]

# Colors
cancer_cols <- c("prostate"="violet", "bladder"="orange")
pt_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> 
  sample(size = n_distinct(meta$patient_deid), replace = F)
names(pt_cols) <- unique(meta$patient_deid)
sub_cols <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]][["Tableau 20"]][["value"]][1:n_distinct(meta$sub_types)]
names(sub_cols) <- unique(meta$sub_types)

# HVGs: Same strategy as before ------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~1)
keep <- rowSums(DESeq2::counts(dds)) >= 100
dds <- dds[keep,]
DESeq2::sizeFactors(dds) <- meta$q_norm_qfactors # Use the Q3 size factors
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- S4Vectors::mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$patient_deid, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = pt_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_patient.pdf", device = "pdf", width = 10, height = 8)

# pdf(file = "heatmap_unsupervised_hclust_1500hvgs_viridis_with_patient_ids.pdf", width = 10, height = 12)
ha <- HeatmapAnnotation(cancer_type = meta$cancer_type, 
                        patient = meta$patient_deid,
                        sub_types = meta$sub_types,
                        col = list("cancer_type"=cancer_cols, "patient"=pt_cols, "sub_types"=sub_cols)
)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.025), quantile(mat, 0.975), length.out=51),
                                   colors = viridis::viridis(51))
)
# dev.off()

# HVGs: regressing out patient effects -----------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~patient_deid)
keep <- rowSums(DESeq2::counts(dds)) >= 100
dds <- dds[keep,]
DESeq2::sizeFactors(dds) <- meta$q_norm_qfactors # Use the Q3 size factors
means <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
dds <- dds[means > 10, ]
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- S4Vectors::mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$patient_deid, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = pt_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_patient_v2.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$cancer_type, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = cancer_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_cancer_type_v2.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$sub_types, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = sub_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_tissue_subtype_v2.pdf", device = "pdf", width = 10, height = 8)

# Heatmap with these 1500 genes
# pdf(file = "heatmap_unsupervised_hclust_1500hvgs_viridis_with_patient_ids_v2.pdf", width = 10, height = 12)
ha <- HeatmapAnnotation(cancer_type = meta$cancer_type, 
                        patient = meta$patient_deid,
                        sub_types = meta$sub_types,
                        col = list("cancer_type"=cancer_cols, "patient"=pt_cols, "sub_types"=sub_cols)
)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.025), quantile(mat, 0.975), length.out=51),
                                   colors = viridis::viridis(51))
)
# dev.off()

# Another interesting strategy ChatGPT helped me come up with ------------------
# 1. Variance-stabilize the data while respecting the design
vsd <- DESeq2::vst(dds, blind = FALSE)  # Keeps patient effect in the transformation

# 2. Extract VST expression matrix
vst_mat <- vsd@assays@data@listData[[1]]

# 3. Model matrix for ~ patient
X <- model.matrix(~patient_deid, data = meta)

# 4. Get residual variance per gene
hvg_resvar <- apply(vst_mat, 1, function(y) {
  fit <- lm.fit(X, y)
  mean(fit$residuals^2)
})

# 5. Select top HVGs
top_hvg_genes <- names(sort(hvg_resvar, decreasing = TRUE))[1:1500]
mat <- vst_mat[top_hvg_genes,] |> apply(1, scale)
# pdf(file = "heatmap_unsupervised_hclust_1500hvgs_viridis_with_patient_ids_v3.pdf", width = 10, height = 12)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nVST Data", 
        show_column_names = F, 
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.025), quantile(mat, 0.975), length.out=51),
                                   colors = viridis::viridis(51))
)
# dev.off()

pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$patient_deid, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = pt_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_patient_v3.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$cancer_type, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = cancer_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_cancer_type_v3.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$sub_types, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = sub_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_tissue_subtype_v3.pdf", device = "pdf", width = 10, height = 8)

# Notes ------------------------------------------------------------------------
# I think that our first strategy should be fine. The other two strategies work 
# well though. The third stragey yields slightly better clustering results in 
# my opinion. It uses VST, not log2(Q3+1). This makes me ask: what is the big 
# whoop with Q3 normalization for GeoMx? It does not seem totally necessary in 
# my opinion. Other strategies can work just as well, if not better.

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
#   [1] stats4    grid      stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] GeomxTools_3.10.0        NanoStringNCTools_1.14.0 S4Vectors_0.44.0         Biobase_2.66.0           BiocGenerics_0.52.0     
# [6] ggbiplot_0.6.2           presto_1.0.0             data.table_1.16.4        Rcpp_1.0.14              lme4_1.1-37             
# [11] Matrix_1.7-3             patchwork_1.3.0          magrittr_2.0.3           ComplexHeatmap_2.22.0    ggplot2_3.5.2           
# [16] dplyr_1.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4                gridExtra_2.3               readxl_1.4.5                rlang_1.1.5                 clue_0.3-66                
# [6] GetoptLong_1.0.5            matrixStats_1.5.0           compiler_4.4.2              reshape2_1.4.4              systemfonts_1.2.2          
# [11] png_0.1-8                   vctrs_0.6.5                 stringr_1.5.1               pkgconfig_2.0.3             shape_1.4.6.1              
# [16] crayon_1.5.3                fastmap_1.2.0               XVector_0.46.0              labeling_0.4.3              UCSC.utils_1.2.0           
# [21] nloptr_2.2.1                ggbeeswarm_0.7.2            purrr_1.0.4                 zlibbioc_1.52.0             GenomeInfoDb_1.42.3        
# [26] jsonlite_2.0.0              EnvStats_3.0.0              DelayedArray_0.32.0         uuid_1.2-1                  BiocParallel_1.40.2        
# [31] parallel_4.4.2              cluster_2.1.6               R6_2.6.1                    stringi_1.8.4               RColorBrewer_1.1-3         
# [36] GGally_2.2.1                limma_3.62.2                parallelly_1.43.0           boot_1.3-31                 numDeriv_2016.8-1.1        
# [41] cellranger_1.1.0            GenomicRanges_1.58.0        scattermore_1.2             SummarizedExperiment_1.36.0 iterators_1.0.14           
# [46] future.apply_1.11.3         IRanges_2.40.1              splines_4.4.2               tidyselect_1.2.1            abind_1.4-8                
# [51] viridis_0.6.5               doParallel_1.0.17           codetools_0.2-20            listenv_0.9.1               lmerTest_3.1-3             
# [56] lattice_0.22-6              tibble_3.2.1                plyr_1.8.9                  withr_3.0.2                 future_1.40.0              
# [61] ggstats_0.9.0               zip_2.3.2                   circlize_0.4.16             Biostrings_2.74.1           pillar_1.10.1              
# [66] BiocManager_1.30.25         MatrixGenerics_1.18.1       renv_1.1.1                  foreach_1.5.2               reformulas_0.4.0           
# [71] generics_0.1.3              sp_2.2-0                    munsell_0.5.1               scales_1.3.0                minqa_1.2.8                
# [76] globals_0.17.0              glue_1.8.0                  pheatmap_1.0.12             tools_4.4.2                 openxlsx_4.2.8             
# [81] locfit_1.5-9.12             ggiraph_0.8.13              dotCall64_1.2               tidyr_1.3.1                 rbibutils_2.3              
# [86] edgeR_4.4.2                 colorspace_2.1-1            nlme_3.1-166                GenomeInfoDbData_1.2.13     beeswarm_0.4.0             
# [91] vipor_0.4.7                 cli_3.6.4                   spam_2.11-1                 ggthemes_5.1.0              S4Arrays_1.6.0             
# [96] viridisLite_0.4.2           gtable_0.3.6                DESeq2_1.46.0               digest_0.6.37               progressr_0.15.1           
# [101] SparseArray_1.6.2           ggrepel_0.9.6               SeuratObject_5.0.2          rjson_0.2.23                htmlwidgets_1.6.4          
# [106] farver_2.1.2                htmltools_0.5.8.1           lifecycle_1.0.4             httr_1.4.7                  GlobalOptions_0.1.2        
# [111] statmod_1.5.0               MASS_7.3-61        

