rm(list = ls())

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/wu-gu-prostate-bladder/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"    

# Packages
library(dplyr)
library(ggplot2)
library(magrittr)
library(consensusMIBC)

# Data 
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta.csv", row.names = 1, )

# Remove NegProbe-WTX before doing any of this analysis
norm <- norm[rownames(norm) != "NegProbe-WTX",]

# Ting was interested in the following classifier designed for Bladder cancer: 
# https://github.com/cit-bioinfo/consensusMIBC

# First, let's subset to only the bladder cancer samples
bladder_norm <- norm[,meta[meta$cancer_type == "bladder",] |> rownames()]

# Trying it out... it seemed to work.
resclass <- consensusMIBC::getConsensusClass(x = bladder_norm, minCor = 0.2, gene_id = "hgnc_symbol")

# I am not sure about the gene symbols. For GeoMx data, some may be unique symbols for certain probes, not standardized symbols.
# Let's check with biomaRt. 
library(biomaRt)
symbols <- rownames(norm)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_key <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", "external_gene_name"), filters = "hgnc_symbol", values = symbols, mart=human)
(symbols %in% gene_key$external_gene_name) |> sum()
# [1] 14830

# There are ~160 genes that have no ensemble accounted for. 

# I tried a few things: 
# 1)
# I printed the missing genes out and passing them to this website: https://biit.cs.ut.ee/gprofiler/convert
# It did not really work.
for (sy in (symbols[!(symbols %in% gene_key$hgnc_symbol)])) {
  cat(sy, "\n")
}

# 2) 
# I searched through the geomx object. GeneID has entrez IDs in it. However, one gene can have multiple IDs.
obj <- readRDS("geomx_obj_processed.RDS")
df <- obj@featureData@data
df |> View()

# 3)
# I exported the initial dataset from DSP control center. It has a key for TargetName to de facto entrez ID.
# But, this does not really work either.
initial <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/wu-gu-prostate-bladder-data/Initial Dataset.xlsx", 
                               sheet = "BioProbeCountMatrix", 
                               cols = c(3, 4, 12))
initial |> group_by(GeneID) |> tally() |> arrange(desc(n))

# norm$enterezid <- plyr::mapvalues(x = rownames(norm), from = initial$TargetName, to = initial$GeneID)
# norm <- relocate(norm, enterezid)
# norm$enterezid %<>% as.numeric()
# rownames(norm) <- NULL
# norm <- tibble::column_to_rownames(norm, var = "enterezid")

# 4) 
# Use the geomx object, and manually look into some genes.
df$entrezid <- df$GeneID |> strsplit(split = ", ") |> sapply(FUN = "[[", 1)
ids <- df |> group_by(entrezid) |> tally() |> arrange(desc(n)) |> filter(n > 1) |> pull(entrezid)

# These are the genes with duplicate entrez ids.
smalldf <- df[df$entrezid %in% ids,]
smalldf |> View()

# We will look into each of them that has multiple possible entrex ids, and choose the one that makes sense.
# Hopefully, over time we will accumulate extrez IDs for all genes in the WTA.
for (entid in smalldf$GeneID) {
  cat(entid, "\n")
}
newentrez <- c(
  2688,
  391194,
  127059,
  401993,
  442186,
  55486,
  51131,
  441457,
  729857,
  100134934,
  200895,
  145553,
  1018,
  442590,
  2993,
  2994,
  5413,
  276,
  5903,
  4738,
  200350,
  3811,
  30014,
  22824,
  30851,
  51390,
  53916,
  494115,
  280,
  1443,
  10214,
  58493,
  64663,
  64682,
  7473,
  81033,
  79541,
  83852,
  89780,
  112398,
  113457,
  728712,
  58508,
  3757,
  5026,
  386676,
  100652748,
  1719,
  7278,
  110599588,
  54529,
  85318,
  26584,
  26581,
  8263,
  349334,
  2812,
  100132285,
  386679,
  25832,
  401427,
  442185,
  343563,
  246721,
  653505,
  728945,
  102723555,
  548313,
  100287932,
  112714,
  474383,
  200030,
  149013,
  101060226,
  101060684,
  54754,
  548644,
  27316,
  100134938,
  107983993,
  100287178,
  100287205
)
smalldf$entrezid <- newentrez
df[df$TargetName %in% smalldf$TargetName,]$entrezid <- smalldf$entrezid

# Now, we actually have unique entrex IDs for all genes.
df |> group_by(entrezid) |> tally() |> arrange(desc(n))

bladder_norm$entrezid <- plyr::mapvalues(x = rownames(bladder_norm), from = df$TargetName, to = df$entrezid)
bladder_norm <- relocate(bladder_norm, entrezid)
rownames(bladder_norm) <- NULL
bladder_norm <- tibble::column_to_rownames(bladder_norm, var = "entrezid")

# Now, let's re-run:
resclassfinal <- consensusMIBC::getConsensusClass(x = bladder_norm, minCor = 0.2, gene_id = "entrezgene")

# It appears that it made a slight difference.
# Now, let's look at the results and what they actually mean. 

# The consensus classification identifies 6 molecular classes: 
#   Luminal Papillary (LumP)
#   Luminal Non Specified (LumNS)
#   Luminal Unstable (LumU)
#   Stroma-rich
#   Basal/Squamous (Ba/Sq)
#   Neuroendocrine-like (NE-like)

# How many of each did we identify? 
resclassfinal$consensusClass |> table()
# Ba/Sq       LumNS        LumP        LumU     NE-like Stroma-rich 
#     9           2          11           7           3           2 

# separationLevel gives a measure of how a sample is representative of its consensus class. 
# It ranges from 0 to 1, with 0 meaning the sample is too close to other consensus classes 
# to be confidently assigned one consensus class label, and 1 meaning the sample is highly 
# representative of its consensus class and well separated from the other consensus classes.

pdf(file = "separation_level_histogram.pdf", width = 6, height = 6)
hist(x = resclassfinal$separationLevel, breaks = 25, main = "separation level", xlab = "value")
dev.off()

pdf(file = "correlation_radar_charts.pdf", width = 12, height = 12)
consensusMIBC::plotCorrelations(xres = resclassfinal)
dev.off()

# What will we do on a case-by-case basis? Do all the samples from the same patient have the same subtype class?
meta$consenusClass <- plyr::mapvalues(x = rownames(meta), from = rownames(resclassfinal), to = resclassfinal$consensusClass)
meta$separationLevel <- plyr::mapvalues(x = rownames(meta), from = rownames(resclassfinal), to = resclassfinal$separationLevel)
bladder_meta <- meta[meta$cancer_type == "bladder",]
table(bladder_meta$patient_deid, bladder_meta$consenusClass)

#       Ba/Sq LumNS LumP LumU NE-like Stroma-rich
# pt_12     0     0    2    0       2           0
# pt_16     1     0    0    0       0           0
# pt_18     0     0    2    1       0           0
# pt_19     2     0    0    0       0           0
# pt_2      0     0    0    0       1           0
# pt_20     0     0    1    0       0           0
# pt_21     0     0    2    0       0           0
# pt_22     0     1    1    0       0           0
# pt_24     0     0    0    1       0           1
# pt_25     0     0    0    1       0           0
# pt_27     0     0    0    0       0           1
# pt_3      1     0    1    0       0           0
# pt_5      0     1    0    1       0           0
# pt_6      0     0    0    1       0           0
# pt_7      0     0    0    2       0           0
# pt_8      4     0    0    0       0           0
# pt_9      1     0    2    0       0           0

# We will save the results and send to Ting.
bladder_meta_simple <- bladder_meta |> 
  dplyr::select(patient_deid, tma_core_number, roi, segment, tags, addn_comments, sub_types, consenusClass, separationLevel) |> 
  arrange(patient_deid, tma_core_number)
classcounts <- table(bladder_meta$patient_deid, bladder_meta$consenusClass) |> as.data.frame.matrix()
classcounts["total",] <- colSums(classcounts)

openxlsx::write.xlsx(x = list("final_result" = resclassfinal, 
                              "bladder_metadata" = bladder_meta_simple, 
                              "class_counts_by_pt" = classcounts), 
                     file = "consensusMIBC_classification_results.xlsx", 
                     rowNames = T)

# Saving the entrez IDs for future use.
write.csv(x = df |> dplyr::select(TargetName, entrezid), file = "TargetName_to_entrezid.csv")

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
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] GeomxTools_3.10.0        NanoStringNCTools_1.14.0 S4Vectors_0.44.0         Biobase_2.66.0           BiocGenerics_0.52.0     
# [6] biomaRt_2.62.1           consensusMIBC_1.1.0      fmsb_0.7.6               magrittr_2.0.3           ggplot2_3.5.2           
# [11] dplyr_1.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            DBI_1.2.3               httr2_1.1.2             readxl_1.4.5            rlang_1.1.5            
# [6] compiler_4.4.2          RSQLite_2.3.9           png_0.1-8               systemfonts_1.2.3       vctrs_0.6.5            
# [11] reshape2_1.4.4          stringr_1.5.1           pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
# [16] dbplyr_2.5.0            XVector_0.46.0          utf8_1.2.4              UCSC.utils_1.2.0        nloptr_2.2.1           
# [21] ggbeeswarm_0.7.2        purrr_1.0.4             bit_4.6.0               zlibbioc_1.52.0         cachem_1.1.0           
# [26] GenomeInfoDb_1.42.3     jsonlite_2.0.0          progress_1.2.3          EnvStats_3.0.0          blob_1.2.4             
# [31] uuid_1.2-1              parallel_4.4.2          prettyunits_1.2.0       R6_2.6.1                stringi_1.8.4          
# [36] RColorBrewer_1.1-3      GGally_2.2.1            parallelly_1.43.0       boot_1.3-31             cellranger_1.1.0       
# [41] numDeriv_2016.8-1.1     Rcpp_1.0.14             future.apply_1.11.3     IRanges_2.40.1          Matrix_1.7-3           
# [46] splines_4.4.2           tidyselect_1.2.1        rstudioapi_0.17.1       codetools_0.2-20        curl_6.2.2             
# [51] listenv_0.9.1           lattice_0.22-6          tibble_3.2.1            lmerTest_3.1-3          plyr_1.8.9             
# [56] withr_3.0.2             KEGGREST_1.46.0         future_1.40.0           ggstats_0.9.0           BiocFileCache_2.14.0   
# [61] zip_2.3.2               xml2_1.3.8              Biostrings_2.74.1       pillar_1.10.1           BiocManager_1.30.25    
# [66] filelock_1.0.3          renv_1.1.1              reformulas_0.4.0        generics_0.1.3          sp_2.2-0               
# [71] hms_1.1.3               munsell_0.5.1           scales_1.3.0            minqa_1.2.8             globals_0.17.0         
# [76] glue_1.8.0              pheatmap_1.0.12         tools_4.4.2             data.table_1.16.4       lme4_1.1-37            
# [81] openxlsx_4.2.8          ggiraph_0.8.13          dotCall64_1.2           grid_4.4.2              tidyr_1.3.1            
# [86] rbibutils_2.3           AnnotationDbi_1.68.0    colorspace_2.1-1        nlme_3.1-166            GenomeInfoDbData_1.2.13
# [91] beeswarm_0.4.0          vipor_0.4.7             cli_3.6.4               rappdirs_0.3.3          spam_2.11-1            
# [96] ggthemes_5.1.0          gtable_0.3.6            digest_0.6.37           progressr_0.15.1        rjson_0.2.23           
# [101] htmlwidgets_1.6.4       SeuratObject_5.0.2      memoise_2.0.1           htmltools_0.5.8.1       lifecycle_1.0.4        
# [106] httr_1.4.7              bit64_4.6.0-1           MASS_7.3-61         

