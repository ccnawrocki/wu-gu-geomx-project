rm(list = ls())

# This script is adapted from the following vignette from NanoString: 
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# All writing and saving functions are commented out to avoid overwriting anything. 

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/wu-gu-prostate-bladder/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Plots
#pdf(file = "QC_plots.pdf", width = 8, height = 8)

# Packages
library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(dplyr)
library(stringr)
library(ggplot2)
library(knitr)
library(reshape2)
library(cowplot)

# GeoMx object
geomx_obj <- readRDS("geomx_obj_unprocessed.RDS")

### QC ###
# QC parameters
ntc_cts <- geomx_obj@protocolData@data[["NTC"]] |> mean() # We know that a couple AOIs have less counts than the NTC well
QC_params <-
  list(minSegmentReads = ntc_cts, # Minimum number of reads (1000) --> We are using the NTC counts as the cutoff
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%) --> default
       percentStitched = 80,   # Minimum % of reads stitched (80%) --> default
       percentAligned = 80,    # Minimum % of reads aligned (80%) --> default
       percentSaturation = 65, # Minimum sequencing saturation (50%) --> Higher than default
       minNegativeCount = 1,   # Minimum negative control counts (10) --> less than default
       maxNTCCount = Inf,     # Maximum counts observed in NTC well (1000) --> We are choosing to ignore this
       minNuclei = 50,         # Minimum # of nuclei estimated (100) --> less than default
       minArea = 500)         # Minimum segment area (500) --> default

# Flagging
geomx_obj <-
  setSegmentQCFlags(geomx_obj, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(geomx_obj)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Graphical summaries of QC statistics
col_by <- "segment"
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = col_by)) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(geomx_obj), "Trimmed (%)", "`slide name`", 80)
QC_histogram(sData(geomx_obj), "Stitched (%)", "`slide name`", 80)
QC_histogram(sData(geomx_obj), "Aligned (%)", "`slide name`", 80)
QC_histogram(sData(geomx_obj), "Saturated (%)", "`slide name`", 65) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(geomx_obj), "area", "`slide name`", 500, scale_trans = "log10")
QC_histogram(sData(geomx_obj), "nuclei", "`slide name`", 50)

# Calculate the negative geometric means for each module
pkcs <- annotation(geomx_obj)
modules <- gsub(".pkc", "", pkcs)
negativeGeoMeans <- 
  esBy(negativeControlSubset(geomx_obj), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(geomx_obj)[["NegGeoMean"]] <- negativeGeoMeans

# Explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(geomx_obj)[, negCols] <- sData(geomx_obj)[["NegGeoMean"]]

# Show the Negative geoMeans plot for each module
for(ann in negCols) {
  plt <- QC_histogram(pData(geomx_obj), ann, "`slide name`", 2, scale_trans = "log10")
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(geomx_obj) <- pData(geomx_obj)[, !colnames(pData(geomx_obj)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(geomx_obj)$NTC),
      col.names = c("NTC Count", "# of Segments"))

# |NTC Count | # of Segments|
# |:---------|-------------:|
# |36923     |            95|

# Manually checking the NTC well... this does not look good, but it is what we have to work with.
ntcwell <- readDccFile(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/wu-gu-prostate-bladder-data/B/dccs/DSP-1001660039570-B-A01.dcc")
ntcwell[["Code_Summary"]][["Count"]] |> hist(breaks = 100)
ntcwell[["Code_Summary"]][["Count"]] |> quantile()
# 0%  25%  50%  75% 100% 
# 1    1    2    3   52 

ntcwell[["Code_Summary"]][["Count"]] |> sum()
# [1] 36923

# Proportion of AOIs that have total probe counts less than the NTC well's counts
mean((geomx_obj@assayData$exprs |> colSums()) < geomx_obj@protocolData@data[["NTC"]])
# [1] 0.01052632

# Show the QC summary... says that 7 AOIs should be filtered, which I think is reasonable.
kable(QC_Summary, caption = "QC Summary Table for each Segment")
# Table: QC Summary Table for each Segment
# 
#   |              | Pass| Warning|
#   |:-------------|----:|-------:|
#   |LowReads      |   94|       1|
#   |LowTrimmed    |   95|       0|
#   |LowStitched   |   95|       0|
#   |LowAligned    |   95|       0|
#   |LowSaturation |   89|       6|
#   |LowNegatives  |   94|       1|
#   |HighNTC       |   95|       0|
#   |LowNuclei     |   94|       1|
#   |LowArea       |   95|       0|
#   |TOTAL FLAGS   |   88|       7|

# Filtering flagged segments out
geomx_obj <- geomx_obj[, QCResults$QCStatus == "PASS"]
dim(geomx_obj)
# Features  Samples 
# 18815       88

# At this point, we are in discussion with NanoString about the NTC well stuff (5/28/2025)

# Removing poor probes
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to
# FALSE if you do not want to remove local outliers
geomx_obj <- setBioProbeQCFlags(geomx_obj,
                                qcCutoffs = list(minProbeRatio = 0.1,
                                                 percentFailGrubbs = 20),
                                removeLocalOutliers = TRUE)

ProbeQCResults <- fData(geomx_obj)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df
# Passed Global Local
# 1  18805      0    10

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <-
  subset(geomx_obj,
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
# Features  Samples 
# 18815       88 
geomx_obj <- ProbeQCPassed

# Check how many unique targets the object has
length(unique(featureData(geomx_obj)[["TargetName"]]))
# [1] 18677

# collapse to gene targets
bygene_geomx_obj <- aggregateCounts(object = geomx_obj)
dim(bygene_geomx_obj)
# Features  Samples 
# 18677       88 

# Looks like I would expect
exprs(bygene_geomx_obj)[1:5, 1:2] 
#       DSP-1001660039570-B-A02.dcc DSP-1001660039570-B-A03.dcc
# A2M                           106                         221
# NAT2                           32                          34
# ACADM                          80                          75
# ACADS                         121                         110
# ACAT1                         182                         157

# Define LOQ SD threshold and minimum value
cutoff <- 1
minLOQ <- 1

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(bygene_geomx_obj))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(bygene_geomx_obj)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(bygene_geomx_obj)[, vars[1]] *
             pData(bygene_geomx_obj)[, vars[2]] ^ cutoff)
  }
}
pData(bygene_geomx_obj)$LOQ <- LOQ

# Filtering based on LOQ
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(bygene_geomx_obj)$Module == module
  Mat_i <- t(esApply(bygene_geomx_obj[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(bygene_geomx_obj)$TargetName, ]

# Save detection rate information to pheno data
pData(bygene_geomx_obj)$GenesDetected <-
  colSums(LOQ_Mat, na.rm = TRUE)
pData(bygene_geomx_obj)$GeneDetectionRate <-
  pData(bygene_geomx_obj)$GenesDetected / nrow(bygene_geomx_obj)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(bygene_geomx_obj)$DetectionThreshold <-
  cut(pData(bygene_geomx_obj)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(bygene_geomx_obj),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = `slide name`)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "`slide name`")
ggplot(pData(bygene_geomx_obj),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segement Type")

# cut percent genes detected at <10%... I do not think anything will get cut here
kable(table(pData(bygene_geomx_obj)$DetectionThreshold,
            pData(bygene_geomx_obj)$segment))
bygene_geomx_obj <-
  bygene_geomx_obj[, pData(bygene_geomx_obj)$GeneDetectionRate >= 0.1]
dim(bygene_geomx_obj)
# Features  Samples 
# 18677       88 

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(bygene_geomx_obj)]
fData(bygene_geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(bygene_geomx_obj)$DetectionRate <-
  fData(bygene_geomx_obj)$DetectedSegments / nrow(pData(bygene_geomx_obj))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(bygene_geomx_obj)[goi, "DetectedSegments"],
  DetectionRate = scales::percent(fData(bygene_geomx_obj)[goi, "DetectionRate"]))
goi_df
# Gene Number DetectionRate
# 1  PDCD1     67         76.1%
# 2  CD274     81         92.0%
# 3   IFNG      3          3.4%
# 4   CD8A     45         51.1%
# 5   CD68     66         75.0%
# 6  EPCAM     73         83.0%
# 7  KRT18     88        100.0% <- this makes sense
# 8  NPHS1      6          6.8%
# 9  NPHS2     22         25.0%
# 10 CALB1     16         18.2%
# 11 CLDN8     60         68.2%

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(bygene_geomx_obj)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(bygene_geomx_obj))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
# Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(bygene_geomx_obj), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
bygene_geomx_obj <-
  bygene_geomx_obj[fData(bygene_geomx_obj)$DetectionRate >= 0.1 |
                     fData(bygene_geomx_obj)$TargetName %in% neg_probes, ]
dim(bygene_geomx_obj)
# Features  Samples 
# 14995       88 

# retain only detected genes of interest
goi <- goi[goi %in% rownames(bygene_geomx_obj)]

# Add the quantile-3-normalization to the object
bygene_geomx_obj <- GeomxTools::normalize(object = bygene_geomx_obj,
                                          norm_method = "quant",
                                          desiredQuantile = .75,
                                          toElt = "q_norm")

### Saving ###
# Saving the counts, the metadata, and the "processed" object
#write.csv(bygene_geomx_obj@assayData$exprs, file = "counts.csv")
meta <- cbind(bygene_geomx_obj@phenoData@data |> dplyr::select(-LOQ),
              sData(bygene_geomx_obj) |> dplyr::select(roi, aoi))
rownames(meta) <- gsub(pattern = "-", replacement = ".", x = rownames(meta))
#write.csv(meta, "meta.csv")
#saveRDS(bygene_geomx_obj, file = "geomx_obj_processed.RDS")

#dev.off()

### Session ###
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
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.3            reshape2_1.4.4           knitr_1.50               magrittr_2.0.3           stringr_1.5.1           
# [6] dplyr_1.1.4              openxlsx_4.2.8           data.table_1.16.4        GeoMxWorkflows_1.12.0    GeomxTools_3.10.0       
# [11] NanoStringNCTools_1.14.0 ggplot2_3.5.2            S4Vectors_0.44.0         Biobase_2.66.0           BiocGenerics_0.52.0     
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            readxl_1.4.5            rlang_1.1.5             compiler_4.4.2          png_0.1-8              
# [6] systemfonts_1.2.2       vctrs_0.6.5             pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
# [11] XVector_0.46.0          labeling_0.4.3          rmarkdown_2.29          UCSC.utils_1.2.0        nloptr_2.2.1           
# [16] ggbeeswarm_0.7.2        purrr_1.0.4             xfun_0.52               zlibbioc_1.52.0         GenomeInfoDb_1.42.3    
# [21] jsonlite_2.0.0          EnvStats_3.0.0          uuid_1.2-1              tweenr_2.0.3            data.tree_1.1.0        
# [26] parallel_4.4.2          R6_2.6.1                stringi_1.8.4           RColorBrewer_1.1-3      reticulate_1.42.0      
# [31] GGally_2.2.1            parallelly_1.43.0       boot_1.3-31             cellranger_1.1.0        numDeriv_2016.8-1.1    
# [36] Rcpp_1.0.14             future.apply_1.11.3     IRanges_2.40.1          igraph_2.1.4            Matrix_1.7-3           
# [41] splines_4.4.2           tidyselect_1.2.1        yaml_2.3.10             codetools_0.2-20        listenv_0.9.1          
# [46] lattice_0.22-6          tibble_3.2.1            lmerTest_3.1-3          plyr_1.8.9              withr_3.0.2            
# [51] askpass_1.2.1           evaluate_1.0.3          Rtsne_0.17              future_1.40.0           ggstats_0.9.0          
# [56] polyclip_1.10-7         zip_2.3.2               Biostrings_2.74.1       pillar_1.10.1           BiocManager_1.30.25    
# [61] renv_1.1.1              reformulas_0.4.0        generics_0.1.3          sp_2.2-0                munsell_0.5.1          
# [66] scales_1.3.0            minqa_1.2.8             BiocStyle_2.34.0        globals_0.17.0          glue_1.8.0             
# [71] pheatmap_1.0.12         tools_4.4.2             lme4_1.1-37             RSpectra_0.16-2         ggiraph_0.8.13         
# [76] dotCall64_1.2           grid_4.4.2              tidyr_1.3.1             rbibutils_2.3           umap_0.2.10.0          
# [81] colorspace_2.1-1        networkD3_0.4.1         nlme_3.1-166            GenomeInfoDbData_1.2.13 beeswarm_0.4.0         
# [86] ggforce_0.4.2           vipor_0.4.7             cli_3.6.4               spam_2.11-1             ggthemes_5.1.0         
# [91] gtable_0.3.6            outliers_0.15           digest_0.6.37           progressr_0.15.1        ggrepel_0.9.6          
# [96] rjson_0.2.23            htmlwidgets_1.6.4       SeuratObject_5.0.2      farver_2.1.2            htmltools_0.5.8.1      
# [101] lifecycle_1.0.4         httr_1.4.7              openssl_2.3.2           MASS_7.3-61     

