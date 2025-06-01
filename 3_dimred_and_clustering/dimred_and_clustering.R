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
meta <- read.csv(file = "meta.csv")

# Remove NegProbe-WTX before doing any of this analysis
norm <- norm[rownames(norm) != "NegProbe-WTX",]
cts <- cts[rownames(cts) != "NegProbe-WTX",]

# Colors
cancer_cols <- c("prostate"="violet", "bladder"="orange")
core_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> 
  sample(size = n_distinct(meta$tma_core_number), replace = F)
names(core_cols) <- unique(meta$tma_core_number)
sub_cols <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][["regular"]][["Tableau 20"]][["value"]][1:n_distinct(meta$sub_types)]
names(sub_cols) <- unique(meta$sub_types)

# Identifying the 1500 with the highest CV
variances <- apply(X = norm, MARGIN = 1, FUN = var)
means <- rowMeans(norm)
cvs <- (variances**0.5)/means
top_genes <- sort(cvs, decreasing = T)[1:1500] |> names()

# PCA
mat <- norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$cancer_type) +
  ggthemes::theme_par() + scale_color_manual(values = cancer_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
ggbiplot(pcaout, var.axes = F, groups = meta$tma_core_number) + 
  ggthemes::theme_par() + scale_color_manual(values = core_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
ggbiplot(pcaout, var.axes = F, groups = meta$sub_types) + 
  ggthemes::theme_par() + scale_color_manual(values = sub_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))

# Heatmap with these 1500 genes
ha <- HeatmapAnnotation(cancer_type = meta$cancer_type, 
                        tma_core = meta$tma_core_number,
                        sub_types = meta$sub_types,
                        col = list("cancer_type"=cancer_cols, "tma_core"=core_cols, "sub_types"=sub_cols)
                        )
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F
)

# Identifying genes with a modeling approach... I think that this is a better strategy so I will do more with this.
(((!(cts%%1 == 0)) |> rowMeans()) != 0) |> which() # Ensuring we only have integer counts
# named integer(0)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~1)
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
sizeFactors(dds) <- meta$q_norm_qfactors # Use the Q3 size factors
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$cancer_type, point.size = 4) +
  ggthemes::theme_par() + scale_color_manual(values = cancer_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_cancer_type.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$tma_core_number, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = core_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_tma_core.pdf", device = "pdf", width = 10, height = 8)
ggbiplot(pcaout, var.axes = F, groups = meta$sub_types, point.size = 4) + 
  ggthemes::theme_par() + scale_color_manual(values = sub_cols) + 
  coord_fixed() + guides(color = guide_legend(override.aes = list(size = 5)))
# ggsave(filename = "PCA_tissue_subtype.pdf", device = "pdf", width = 10, height = 8)

# Heatmap with these 1500 genes
ha <- HeatmapAnnotation(cancer_type = meta$cancer_type, 
                        tma_core = meta$tma_core_number,
                        sub_types = meta$sub_types,
                        col = list("cancer_type"=cancer_cols, "tma_core"=core_cols, "sub_types"=sub_cols)
)
pdf(file = "heatmap_unsupervised_hclust_1500hvgs.pdf", width = 10, height = 12)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F, 
)
dev.off()

# Unsupervised, but split into 2 groups
pdf(file = "heatmap_unsupervised_hclust_1500hvgs_2clusters.pdf", width = 10, height = 12)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_row_names = F, 
        show_column_names = F, 
        column_split = 2
)
dev.off()

# Unsupervised, but split into 3 groups
pdf(file = "heatmap_unsupervised_hclust_1500hvgs_3clusters.pdf", width = 10, height = 12)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_row_names = F, 
        show_column_names = F, 
        column_split = 3
)
dev.off()

# Supervised
pdf(file = "heatmap_supervised_hclust_1500hvgs.pdf", width = 10, height = 12)
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_row_names = F, 
        column_split = meta$cancer_type
)
dev.off()

# Notes: 
# It looks pretty good. When we get the patient IDs, we can re-cluster while
# accounting for that as a batch variable for good measure.


# Session
sessionInfo()

# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.4.1
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
#   [1] grid      stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.0          ggbiplot_0.6.2           ComplexHeatmap_2.22.0    cowplot_1.1.3            reshape2_1.4.4          
# [6] knitr_1.50               magrittr_2.0.3           stringr_1.5.1            dplyr_1.1.4              openxlsx_4.2.8          
# [11] data.table_1.16.4        GeoMxWorkflows_1.12.0    GeomxTools_3.10.0        NanoStringNCTools_1.14.0 ggplot2_3.5.2           
# [16] S4Vectors_0.44.0         Biobase_2.66.0           BiocGenerics_0.52.0     
# 
# loaded via a namespace (and not attached):
#   [1] matrixStats_1.5.0           spatstat.sparse_3.1-0       httr_1.4.7                  RColorBrewer_1.1-3         
# [5] doParallel_1.0.17           numDeriv_2016.8-1.1         tools_4.4.2                 sctransform_0.4.1          
# [9] utf8_1.2.4                  R6_2.6.1                    lazyeval_0.2.2              uwot_0.2.3                 
# [13] GetoptLong_1.0.5            withr_3.0.2                 sp_2.2-0                    GGally_2.2.1               
# [17] gridExtra_2.3               progressr_0.15.1            cli_3.6.4                   spatstat.explore_3.4-2     
# [21] fastDummies_1.7.5           labeling_0.4.3              Seurat_5.2.1                spatstat.data_3.1-6        
# [25] ggridges_0.5.6              pbapply_1.7-2               askpass_1.2.1               systemfonts_1.2.2          
# [29] parallelly_1.43.0           readxl_1.4.5                generics_0.1.3              shape_1.4.6.1              
# [33] ica_1.0-3                   spatstat.random_3.3-3       zip_2.3.2                   Matrix_1.7-3               
# [37] ggbeeswarm_0.7.2            abind_1.4-8                 lifecycle_1.0.4             yaml_2.3.10                
# [41] SummarizedExperiment_1.36.0 SparseArray_1.6.2           Rtsne_0.17                  glmGamPoi_1.18.0           
# [45] promises_1.3.2              crayon_1.5.3                miniUI_0.1.2                lattice_0.22-6             
# [49] pillar_1.10.1               GenomicRanges_1.58.0        rjson_0.2.23                boot_1.3-31                
# [53] future.apply_1.11.3         codetools_0.2-20            glue_1.8.0                  ggiraph_0.8.13             
# [57] outliers_0.15               spatstat.univar_3.1-2       vctrs_0.6.5                 png_0.1-8                  
# [61] spam_2.11-1                 Rdpack_2.6.4                cellranger_1.1.0            gtable_0.3.6               
# [65] xfun_0.52                   rbibutils_2.3               S4Arrays_1.6.0              mime_0.13                  
# [69] reformulas_0.4.0            survival_3.7-0              pheatmap_1.0.12             iterators_1.0.14           
# [73] fitdistrplus_1.2-2          ROCR_1.0-11                 nlme_3.1-166                EnvStats_3.0.0             
# [77] RcppAnnoy_0.0.22            GenomeInfoDb_1.42.3         data.tree_1.1.0             irlba_2.3.5.1              
# [81] vipor_0.4.7                 KernSmooth_2.23-24          colorspace_2.1-1            DESeq2_1.46.0              
# [85] tidyselect_1.2.1            compiler_4.4.2              DelayedArray_0.32.0         plotly_4.10.4              
# [89] scales_1.3.0                lmtest_0.9-40               digest_0.6.37               goftest_1.2-3              
# [93] spatstat.utils_3.1-3        minqa_1.2.8                 rmarkdown_2.29              XVector_0.46.0             
# [97] htmltools_0.5.8.1           pkgconfig_2.0.3             lme4_1.1-37                 umap_0.2.10.0              
# [101] MatrixGenerics_1.18.1       fastmap_1.2.0               rlang_1.1.5                 GlobalOptions_0.1.2        
# [105] htmlwidgets_1.6.4           ggthemes_5.1.0              UCSC.utils_1.2.0            shiny_1.10.0               
# [109] farver_2.1.2                zoo_1.8-14                  jsonlite_2.0.0              BiocParallel_1.40.2        
# [113] GenomeInfoDbData_1.2.13     dotCall64_1.2               munsell_0.5.1               Rcpp_1.0.14                
# [117] viridis_0.6.5               reticulate_1.42.0           stringi_1.8.4               zlibbioc_1.52.0            
# [121] MASS_7.3-61                 plyr_1.8.9                  ggstats_0.9.0               parallel_4.4.2             
# [125] listenv_0.9.1               ggrepel_0.9.6               deldir_2.0-4                Biostrings_2.74.1          
# [129] splines_4.4.2               tensor_1.5                  circlize_0.4.16             locfit_1.5-9.12            
# [133] igraph_2.1.4                uuid_1.2-1                  spatstat.geom_3.3-6         RcppHNSW_0.6.0             
# [137] evaluate_1.0.3              SeuratObject_5.0.2          renv_1.1.1                  BiocManager_1.30.25        
# [141] nloptr_2.2.1                foreach_1.5.2               tweenr_2.0.3                httpuv_1.6.16              
# [145] networkD3_0.4.1             RANN_2.6.2                  tidyr_1.3.1                 openssl_2.3.2              
# [149] purrr_1.0.4                 polyclip_1.10-7             future_1.40.0               clue_0.3-66                
# [153] scattermore_1.2             ggforce_0.4.2               xtable_1.8-4                RSpectra_0.16-2            
# [157] later_1.4.2                 viridisLite_0.4.2           tibble_3.2.1                lmerTest_3.1-3             
# [161] beeswarm_0.4.0              IRanges_2.40.1              cluster_2.1.6               globals_0.17.0             
# [165] BiocStyle_2.34.0           

