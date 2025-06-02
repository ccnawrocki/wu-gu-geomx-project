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

# Data 
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta.csv")

# Remove NegProbe-WTX before doing any of this analysis
norm <- norm[rownames(norm) != "NegProbe-WTX",]
cts <- cts[rownames(cts) != "NegProbe-WTX",]

# Need to pseudobulk or use a mixed model to account for intrapatient correlation/pseudo-replication bias. 
# For this, we need the patient identifier for each AOI.

# Adding patient IDs to the data


# Saving the updated metadata


# We will try the mixed model approach first, using a random intercept model with lme4. 
# If we get singular fits, then we will pseudo-bulk by patient and use DESeq2. 

# First, we need a function to do the work. 














