# Packages
library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(data.table)
library(openxlsx)
library(dplyr)
library(stringr)
library(magrittr)

# Defining where the data is
# Using the second sequencing library, denoted as "B"
datadir <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/wu-gu-prostate-bladder-data/B" 

# The directory should have the lab worksheets and the GeoMx_NGS_Pipeline_DCC folders files in it
list.files(datadir) 
# [1] "GeoMx_NGS_Pipeline_DCC"                     "Hs_R_NGS_WTA_v1.0.pkc"                     
# [3] "Wu_Gu Run 1_20250502T1756_LabWorksheet.txt"

# Create the metadata excel file. 
labworksheets <- list.files(datadir)[list.files(datadir) |> grep(pattern = "LabWorksheet")]
metas <- list()
for (sheet in labworksheets) {
  
  # Getting the run name via the worksheet
  run_name <- str_split(string = sheet, pattern = "_[0-9]{3}", simplify = T)[,1]
  
  # Extracting the header and the NTC line from the worksheet
  header <- readLines(con = file.path(datadir, sheet), n = 17) 
  
  # Defining the column names for the metadata
  column_names <- header[16] |> 
    str_split(pattern = "\t") |> 
    unlist()
  
  # Reading the metadata, without the header and the NTC line
  meta <- fread(file.path(datadir, sheet), 
                skip = 17, sep = "\t", header = F) |> 
    as.data.frame()
  
  # Adding the column names
  colnames(meta) <- column_names
  
  # Reformatting the Roi column
  meta$Roi <- gsub(x = meta$Roi, pattern = "=|\"", replacement = "") 
  
  # Adding the NTC line at the top of the metadata
  ntcline <- header[17] |> str_split(pattern = "\t") |> unlist()
  
  # It needs to be the same number of columns as the other lines
  ntcline <- c(ntcline, rep("", length(column_names)-2)) |>
    matrix(nrow = 1, byrow = T, dimnames = list("0", colnames(meta)))
  
  # It needs to have the panel information added to it
  ntcline[which(colnames(meta) == "Panel")] <- meta$Panel[1]
  
  # Actually combining the NTC line and the other metadata
  meta <- rbind(ntcline, meta)
  
  # Adding to the metadata list
  metas[[run_name]] <- meta
}

# Combining all the metadata data frames
metadata <- bind_rows(metas)

# The columns must be in this order
metadata %<>% dplyr::select(Sample_ID, Panel, `Slide Name`, 
                            Roi, Segment, Aoi, Area, Nuclei, 
                            `Scan Width`, `Scan Height`, 
                            `ROI Coordinate X`, `ROI Coordinate Y`, 
                            `Scan Offset X`, `Scan Offset Y`, Tags)

# Some columns must be renamed
metadata %<>% rename(panel=Panel, `slide name`=`Slide Name`, roi=Roi, aoi=Aoi, segment=Segment, nuclei=Nuclei, area=Area)


write.xlsx(metadata, file = "metadata.xlsx")

# Move all the dcc files into one "dcc" folder. Have to do this in bash.
# $ cd ~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/wu-gu-prostate-bladder-data/B
# $ mkdir dccs
# $ cp GeoMx_NGS_Pipeline_DCC*/* dccs

# Creating the GeoMx objects
PKCFiles <- list.files(datadir)[list.files(datadir) |> grep(pattern = ".pkc")]
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
ANNOFile <- "metadata.xlsx"
geomx_obj <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                    pkcFiles = file.path(datadir, PKCFiles),
                                    phenoDataFile = ANNOFile,
                                    phenoDataSheet = "Sheet 1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("aoi", "roi"),
                                    experimentDataColNames = c("panel")
)

# Checking the data looks right
geomx_obj@assayData$exprs[1:5,1:5]
geomx_obj@phenoData@data |> dplyr::glimpse()

# Exporting the geomx object as an RDS file
saveRDS(object = geomx_obj, file = "geomx_obj_unprocessed.RDS")

# Notes: 
# The NTC well had high counts in library A. So, library B was generated. Still, there are high NTC well counts. 
# We are sticking with B, since the NTC well counts are lower. 
# It is likely not a huge problem for this study, since there is only 1 plate. Thus, everything is probably 
# contaminated by the same things/to the same extent. 
# Could we combine libraries A and B to get better sequencing depth? 

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
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] magrittr_2.0.3           stringr_1.5.1            dplyr_1.1.4              openxlsx_4.2.8           data.table_1.16.4       
# [6] GeoMxWorkflows_1.12.0    GeomxTools_3.10.0        NanoStringNCTools_1.14.0 ggplot2_3.5.2            S4Vectors_0.44.0        
# [11] Biobase_2.66.0           BiocGenerics_0.52.0     
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.4            readxl_1.4.5            rlang_1.1.5             compiler_4.4.2          png_0.1-8              
# [6] systemfonts_1.2.2       vctrs_0.6.5             reshape2_1.4.4          pkgconfig_2.0.3         crayon_1.5.3           
# [11] fastmap_1.2.0           XVector_0.46.0          rmarkdown_2.29          UCSC.utils_1.2.0        nloptr_2.2.1           
# [16] ggbeeswarm_0.7.2        purrr_1.0.4             xfun_0.52               zlibbioc_1.52.0         GenomeInfoDb_1.42.3    
# [21] jsonlite_2.0.0          EnvStats_3.0.0          uuid_1.2-1              tweenr_2.0.3            data.tree_1.1.0        
# [26] parallel_4.4.2          R6_2.6.1                stringi_1.8.4           RColorBrewer_1.1-3      reticulate_1.42.0      
# [31] GGally_2.2.1            parallelly_1.43.0       boot_1.3-31             cellranger_1.1.0        numDeriv_2016.8-1.1    
# [36] knitr_1.50              Rcpp_1.0.14             future.apply_1.11.3     IRanges_2.40.1          igraph_2.1.4           
# [41] Matrix_1.7-3            splines_4.4.2           tidyselect_1.2.1        yaml_2.3.10             codetools_0.2-20       
# [46] listenv_0.9.1           lattice_0.22-6          tibble_3.2.1            lmerTest_3.1-3          plyr_1.8.9             
# [51] withr_3.0.2             askpass_1.2.1           evaluate_1.0.3          Rtsne_0.17              future_1.40.0          
# [56] ggstats_0.9.0           polyclip_1.10-7         zip_2.3.2               Biostrings_2.74.1       pillar_1.10.1          
# [61] BiocManager_1.30.25     renv_1.1.1              reformulas_0.4.0        generics_0.1.3          sp_2.2-0               
# [66] munsell_0.5.1           scales_1.3.0            minqa_1.2.8             BiocStyle_2.34.0        globals_0.17.0         
# [71] glue_1.8.0              pheatmap_1.0.12         tools_4.4.2             lme4_1.1-37             RSpectra_0.16-2        
# [76] ggiraph_0.8.13          dotCall64_1.2           cowplot_1.1.3           grid_4.4.2              tidyr_1.3.1            
# [81] rbibutils_2.3           umap_0.2.10.0           colorspace_2.1-1        networkD3_0.4.1         nlme_3.1-166           
# [86] GenomeInfoDbData_1.2.13 beeswarm_0.4.0          ggforce_0.4.2           vipor_0.4.7             cli_3.6.4              
# [91] spam_2.11-1             ggthemes_5.1.0          gtable_0.3.6            outliers_0.15           digest_0.6.37          
# [96] progressr_0.15.1        ggrepel_0.9.6           rjson_0.2.23            htmlwidgets_1.6.4       SeuratObject_5.0.2     
# [101] farver_2.1.2            htmltools_0.5.8.1       lifecycle_1.0.4         httr_1.4.7              openssl_2.3.2          
# [106] MASS_7.3-61

