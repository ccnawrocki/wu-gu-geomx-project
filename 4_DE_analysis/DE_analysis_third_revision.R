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

# Data 
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta_amended.csv", row.names = 1)
norm <- norm[rownames(norm) != "NegProbe-WTX", rownames(meta)]
cts <- cts[rownames(cts) != "NegProbe-WTX", rownames(meta)]

# This analysis only looks at the prostate data
prostate <- meta[(meta$cancer_type == "prostate"),] |> rownames()
cts <- cts[,prostate]
norm <- norm[,prostate]
meta <- meta[prostate,]

# For visualization
snorm <- apply(X = norm, MARGIN = 1, FUN = scale) |> t()
colnames(snorm) <- sprintf("%03d", meta$roi)
bc <- limma::removeBatchEffect(x = norm, batch = meta$patient_deid, design = model.matrix(~sub_types_v2, data = meta))
sbc <- apply(X = bc, MARGIN = 1, FUN = scale) |> t()
colnames(sbc) <- sprintf("%03d", meta$roi)

# Q11: could you please delete “Acinar IDC-P groups (2 AOIs)” as in other comparisons, and re-run the plot diagram and heat map below? Thank you!
meta$sub_types_v5 <- case_when(meta$roi %in% c(48, 49) ~ "Precursor", 
                               meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "AIP") ~ "Intraductal_spread", 
                               T ~ meta$sub_types_v2)

small_meta <- meta |> filter(sub_types_v5 %in% c("Precursor", "Intraductal_spread"))
small_meta$group <- small_meta$sub_types_v5

# Calculating the difference in means across the two groups for every gene.
mm <- model.matrix(~0+group, data = small_meta)
small_norm <- norm[,rownames(small_meta)]
groupmeans <- (as.matrix(small_norm) %*% mm) |> sweep(x = _, MARGIN = 2, STATS = colSums(mm), FUN = "/")
colnames(groupmeans) <- gsub(pattern = "group", replacement = "", x = colnames(groupmeans))
groupmeans <- as.data.frame(groupmeans)

toplot <- small_norm[which((abs(groupmeans[,"Precursor"]-groupmeans[,"Intraductal_spread"])) > 2),] |> rownames()
idx <- sprintf("%03d", small_meta$roi)
mat <- snorm[toplot, idx]
ha <- HeatmapAnnotation(group = small_meta$group,
                        sub_type = small_meta$sub_types_v2,
                        col = list("group"=c("Precursor"="dodgerblue", "Intraductal_spread"="gold3"), 
                                   "sub_type"=c("Acinar IDC-P_crib"="blue", "AIP"="purple")))
pdf(file = "2.3 - Precursor_vs_Intraductal_spread_top_genes_by_logFC_heatmap.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = small_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "row−scaled\nlog2(q3norm+1)",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

mat <- sbc[toplot, idx]
pdf(file = "2.3 - Precursor_vs_Intraductal_spread_top_genes_by_logFC_heatmap_batch_corrected.pdf", width = 8, height = 14)
Heatmap(matrix = mat, column_split = small_meta$group, top_annotation = ha,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7),
        name = "Batch\ncorrected\nexpression",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(1.5, "mm"), 
        column_title_gp = gpar(col = NA)
)
dev.off()

groupmeans$logFC <- (groupmeans[,"Precursor"] - groupmeans[,"Intraductal_spread"])
write.csv(x = groupmeans, file = "2.3 - Precursor_vs_Intraductal_spread_top_genes_by_logFC.csv")

