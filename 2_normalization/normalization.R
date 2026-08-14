rm(list = ls())

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/wu-gu-prostate-bladder/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(dplyr)
library(ggplot2)
library(magrittr)
library(patchwork)

### Data ### -------------------------------------------------------------------

# Counts and metadata
cts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("meta.csv", row.names = 1)
all(rownames(meta) == colnames(cts))
# [1] TRUE

# Adding to metadata 
samplesheet <- openxlsx::read.xlsx("Wu_Gu_Run1_LabWorksheet.xlsx")
samplesheet$Sample_ID <- gsub(pattern = "\\-A\\-", replacement = "\\-B\\-", x = samplesheet$Sample_ID)
samplesheet$Sample_ID <- gsub(pattern = "\\-", replacement = "\\.", x = samplesheet$Sample_ID)
samplesheet$Sample_ID <- paste(samplesheet$Sample_ID, "dcc", sep = ".")
meta$Sample_ID <- rownames(meta)
meta <- inner_join(x = meta, 
                   y = samplesheet |> select(Sample_ID, TMA.core.number, `Sub-types`, Addn.comments), 
                   by = "Sample_ID")

# Organizing metadata a bit more
rownames(meta) <- meta$Sample_ID
colnames(meta) %<>% tolower()
colnames(meta) <- gsub(pattern = "\\.|\\-", replacement = "_", x = colnames(meta))
meta$cancer_type <- ifelse(test = meta$slide_name == "Wu/Gu TMA 1", yes = "prostate", no = "bladder")

# Re-saving the metadata
# write.csv(x = meta, file = "meta.csv")

### Q3 Normalization ### -------------------------------------------------------
# This is equivalent to the Q3 normalization in the GeomxTools package. 
qs <- apply(X = cts, MARGIN = 2, FUN = quantile, 0.75)
nfs <- qs / EnvStats::geoMean(qs)
q3norm <- sweep(x = cts, MARGIN = 2, STATS = nfs, FUN = "/")

# Saving Q3 normalized data
# write.csv(x = q3norm, file = "q3norm.csv")

# I will be working with the log-transformation of the Q3 normalized data.
lq3 <- log2(q3norm + 1)

# Saving the logQ3 data as well
# write.csv(x = lq3, file = "log2plus1_q3norm.csv")

# Computing counts per million (CPM) and log2(CPM+1) as well for comparison
cpm <- sweep(x = cts, MARGIN = 2, STATS = colSums(cts), FUN = "/") * 1e6
lcpm <- log2(cpm + 1)

# Viz
# pdf(file = "normalization_plots.pdf", width = 10, height = 10)
boxplot(cts+1,
        log = "y",
        col = "red3", 
        main = "Raw Counts + 1", 
        names = 1:88, xlab = "AOI #")
boxplot(q3norm+1,
        log = "y",
        col = "dodgerblue", 
        main = "Q3 Normalized Counts + 1", 
        names = 1:88, xlab = "AOI #")
boxplot(lq3,
        col = "gold", 
        main = "log2(Q3+1) Transformed Counts", 
        names = 1:88, xlab = "AOI #")
boxplot(cpm+1,
        log = "y",
        col = "limegreen", 
        main = "CPM Normalized Counts + 1", 
        names = 1:88, xlab = "AOI #")
boxplot(lcpm,
        col = "purple", 
        main = "log2(CPM+1) Normalized Counts", 
        names = 1:88, xlab = "AOI #")
# dev.off()

# Base plotting makes plots that are sort of big. Making some smaller ggplots too:
# pdf(file = "normalization_plots_smaller.pdf", width = 10, height = 10)
ggplot() + 
  geom_boxplot(data = cts |> tidyr::pivot_longer(cols = 1:ncol(cts), names_to = "aoi", values_to = "RawCounts"), 
               mapping = aes(x = aoi, y = RawCounts+1), fill = "red3", outliers = F) + 
  scale_y_log10() + 
  ggthemes::theme_par() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  labs(title = "Raw Counts + 1")
ggplot() + 
  geom_boxplot(data = q3norm |> tidyr::pivot_longer(cols = 1:ncol(q3norm), names_to = "aoi", values_to = "Q3"), 
               mapping = aes(x = aoi, y = Q3+1), fill = "dodgerblue", outliers = F) + 
  scale_y_log10() + 
  ggthemes::theme_par() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  labs(title = "Q3 Normalized Counts + 1")
ggplot() + 
  geom_boxplot(data = lq3 |> tidyr::pivot_longer(cols = 1:ncol(lq3), names_to = "aoi", values_to = "log2(Q3+1)"), 
               mapping = aes(x = aoi, y = `log2(Q3+1)`), fill = "gold", outliers = F) + 
  ggthemes::theme_par() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  labs(title = "log2(Q3+1) Transformed Counts")
ggplot() + 
  geom_boxplot(data = cpm |> tidyr::pivot_longer(cols = 1:ncol(cpm), names_to = "aoi", values_to = "cpm"), 
               mapping = aes(x = aoi, y = cpm+1), fill = "limegreen", outliers = F) + 
  scale_y_log10() + 
  ggthemes::theme_par() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  labs(title = "CPM+1")
ggplot() + 
  geom_boxplot(data = lcpm |> tidyr::pivot_longer(cols = 1:ncol(lcpm), names_to = "aoi", values_to = "log2(CPM+1)"), 
               mapping = aes(x = aoi, y = `log2(CPM+1)`), fill = "purple", outliers = F) + 
  ggthemes::theme_par() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + 
  labs(title = "log2(CPM+1) Transformed Counts")
# dev.off()

### Session ### ----------------------------------------------------------------
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
#   [1] RcppAnnoy_0.0.22        splines_4.4.2           later_1.4.2             tibble_3.2.1            cellranger_1.1.0       
# [6] polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.4         Rdpack_2.6.4            doParallel_1.0.17      
# [11] globals_0.17.0          lattice_0.22-6          MASS_7.3-61             plotly_4.10.4           rmarkdown_2.29         
# [16] yaml_2.3.10             httpuv_1.6.16           Seurat_5.2.1            sctransform_0.4.1       spam_2.11-1            
# [21] zip_2.3.2               askpass_1.2.1           spatstat.sparse_3.1-0   sp_2.2-0                reticulate_1.42.0      
# [26] pbapply_1.7-2           minqa_1.2.8             RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.52.0        
# [31] EnvStats_3.0.0          Rtsne_0.17              purrr_1.0.4             tweenr_2.0.3            circlize_0.4.16        
# [36] GenomeInfoDbData_1.2.13 IRanges_2.40.1          data.tree_1.1.0         ggrepel_0.9.6           irlba_2.3.5.1          
# [41] spatstat.utils_3.1-3    listenv_0.9.1           pheatmap_1.0.12         BiocStyle_2.34.0        umap_0.2.10.0          
# [46] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-3   fitdistrplus_1.2-2      parallelly_1.43.0      
# [51] codetools_0.2-20        ggforce_0.4.2           tidyselect_1.2.1        shape_1.4.6.1           outliers_0.15          
# [56] UCSC.utils_1.2.0        farver_2.1.2            lme4_1.1-37             spatstat.explore_3.4-2  matrixStats_1.5.0      
# [61] jsonlite_2.0.0          GetoptLong_1.0.5        progressr_0.15.1        ggridges_0.5.6          survival_3.7-0         
# [66] iterators_1.0.14        systemfonts_1.2.2       foreach_1.5.2           tools_4.4.2             ica_1.0-3              
# [71] Rcpp_1.0.14             glue_1.8.0              gridExtra_2.3           xfun_0.52               ggthemes_5.1.0         
# [76] GenomeInfoDb_1.42.3     withr_3.0.2             numDeriv_2016.8-1.1     BiocManager_1.30.25     fastmap_1.2.0          
# [81] GGally_2.2.1            boot_1.3-31             openssl_2.3.2           digest_0.6.37           mime_0.13              
# [86] R6_2.6.1                colorspace_2.1-1        networkD3_0.4.1         scattermore_1.2         tensor_1.5             
# [91] spatstat.data_3.1-6     tidyr_1.3.1             generics_0.1.3          renv_1.1.1              httr_1.4.7             
# [96] htmlwidgets_1.6.4       ggstats_0.9.0           uwot_0.2.3              pkgconfig_2.0.3         gtable_0.3.6           
# [101] lmtest_0.9-40           XVector_0.46.0          htmltools_0.5.8.1       dotCall64_1.2           clue_0.3-66            
# [106] SeuratObject_5.0.2      scales_1.3.0            png_0.1-8               spatstat.univar_3.1-2   reformulas_0.4.0       
# [111] rjson_0.2.23            uuid_1.2-1              nlme_3.1-166            nloptr_2.2.1            zoo_1.8-14             
# [116] GlobalOptions_0.1.2     KernSmooth_2.23-24      parallel_4.4.2          miniUI_0.1.2            vipor_0.4.7            
# [121] pillar_1.10.1           vctrs_0.6.5             RANN_2.6.2              promises_1.3.2          xtable_1.8-4           
# [126] cluster_2.1.6           beeswarm_0.4.0          evaluate_1.0.3          cli_3.6.4               compiler_4.4.2         
# [131] rlang_1.1.5             crayon_1.5.3            future.apply_1.11.3     labeling_0.4.3          plyr_1.8.9             
# [136] ggbeeswarm_0.7.2        ggiraph_0.8.13          stringi_1.8.4           deldir_2.0-4            viridisLite_0.4.2      
# [141] lmerTest_3.1-3          munsell_0.5.1           Biostrings_2.74.1       lazyeval_0.2.2          spatstat.geom_3.3-6    
# [146] Matrix_1.7-3            RcppHNSW_0.6.0          future_1.40.0           shiny_1.10.0            rbibutils_2.3          
# [151] ROCR_1.0-11             igraph_2.1.4            readxl_1.4.5 

