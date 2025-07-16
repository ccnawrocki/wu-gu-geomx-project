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

## PROSTATE (TMA1) -------------------------------------------------------------
prostate_meta <- meta |> filter(cancer_type == "prostate")
prostate_cts <- cts[,rownames(prostate_meta)]
prostate_norm <- norm[,rownames(prostate_meta)] 

# Colors
subtype_cols <- InSituType::colorCellTypes(names = unique(prostate_meta$sub_types), palette = "tableau20")
patient_cols <- ggprism::ggprism_data$colour_palettes$colors
names(patient_cols) <- unique(prostate_meta$patient_deid)

# HVGs, using DESeq2 to estimate dispersion# HVGs, usiprism_light2ng DESeq2 to estimate dispersion
dds <- DESeq2::DESeqDataSetFromMatrix(countData = prostate_cts[rownames(prostate_cts) != "NegProbe-WTX",], 
                                      colData = prostate_meta, design = ~1)
DESeq2::sizeFactors(dds) <- prostate_meta$q_norm_qfactors # Use the Q3 size factors
keep <- rowMeans(DESeq2::counts(dds, normalized = T)) >= 10
dds <- dds[keep,]
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- S4Vectors::mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- prostate_norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_patient.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_sub_type.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:15], n_neighbors = 5)
prostate_meta[,c("umap_1", "umap_2")] <- um
ggplot() + 
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = patient_cols) + 
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_patient.pdf", height = 6, width = 8)
ggplot() + 
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_sub_type.pdf", height = 6, width = 8)

# Clearly, the patient identity drives the clustering. Therefore, I will try
# to perform batch-correction with a few techniques.
# This is based on the following sources:
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-after-vst-are-there-still-batches-in-the-pca-plot
# https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html#batch-correction
# https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html#ruvr-estimating-the-factors-of-unwanted-variation-using-residuals

### LIMMA WORKED BEST ###

# ## ComBat_Seq
# combat_bc <- sva::ComBat_seq(counts = prostate_cts[rownames(prostate_cts != "NegProbe-WTX"),] |> as.matrix(), 
#                              batch = prostate_meta$patient_deid, group = prostate_meta$sub_types)
# qs <- apply(X = combat_bc, MARGIN = 2, FUN = quantile, 0.75)
# nfs <- qs / EnvStats::geoMean(qs)
# combat_bc_norm <- log2(sweep(x = combat_bc, MARGIN = 2, STATS = nfs, FUN = "/") + 1)
# 
# ## RUVr
# mm <- model.matrix(~sub_types, data=prostate_meta)
# y <- edgeR::DGEList(counts = prostate_cts[rownames(prostate_cts != "NegProbe-WTX"),] |> as.matrix(), 
#                     samples = prostate_meta)
# y$samples$norm.factors <- y$samples$q_norm_qfactors
# y <- edgeR::estimateGLMCommonDisp(y, mm)
# y <- edgeR::estimateGLMTagwiseDisp(y, mm)
# fit <- edgeR::glmFit(y, mm)
# res <- residuals(fit, type = "deviance")
# ruvr_out <- RUVSeq::RUVr(x = as.matrix(prostate_cts), rownames(prostate_cts), k = 3, res)
# qs <- apply(X = ruvr_out$normalizedCounts, MARGIN = 2, FUN = quantile, 0.75)
# nfs <- qs / EnvStats::geoMean(qs)
# ruvr_bc_norm <- log2(sweep(x = ruvr_out$normalizedCounts, MARGIN = 2, STATS = nfs, FUN = "/") + 1)
# 
# ## RUVg
# ruvg_out <- RUVg(x = as.matrix(prostate_cts) |> round(), cIdx = "NegProbe-WTX", k = 3, isLog = F)
# qs <- apply(X = ruvg_out$normalizedCounts, MARGIN = 2, FUN = quantile, 0.75)
# nfs <- qs / EnvStats::geoMean(qs)
# ruvg_bc_norm <- log2(sweep(x = ruvg_out$normalizedCounts, MARGIN = 2, STATS = nfs, FUN = "/") + 1)

## limma
mm <- model.matrix(~sub_types, data=prostate_meta)
limma_bc <- limma::removeBatchEffect(x = prostate_norm, batch = prostate_meta$patient_deid, design = mm)

## Finding 1500 variable genes, using non-batch-corrected data, modeling with patient, and sorting by MSE
X <- model.matrix(~patient_deid, data = prostate_meta)
hvg_resvar <- apply(prostate_norm, 1, function(y) {
  fit <- lm.fit(X, y)
  mean(fit$residuals^2)
})
top_genes <- names(sort(hvg_resvar, decreasing = TRUE))[1:1500]

# ## Re-doing PCA and UMAP with ComBat_Seq results
# mat <- combat_bc_norm[top_genes,] |> apply(1, scale)
# 
# pcaout <- prcomp(mat, center = F, scale. = F)
# plot(pcaout)
# pca_var<-pcaout$sdev**2
# pve <- pca_var/sum(pca_var)
# plot(pve)
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# um <- uwot::umap(X = pcaout$x[,1:8], n_neighbors = 5, metric = "cosine", min_dist = 0.01, spread = 5)
# prostate_meta[,c("umap_1", "umap_2")] <- um
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))
# 
# ## Re-doing PCA and UMAP with RUVg results
# mat <- ruvg_bc_norm[top_genes,] |> apply(1, scale)
# 
# pcaout <- prcomp(mat, center = F, scale. = F)
# plot(pcaout)
# pca_var<-pcaout$sdev**2
# pve <- pca_var/sum(pca_var)
# plot(pve)
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# um <- uwot::umap(X = pcaout$x[,1:10], n_neighbors = 5, metric = "cosine", min_dist = 0.01, spread = 5)
# prostate_meta[,c("umap_1", "umap_2")] <- um
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))
# 
# ## Re-doing PCA and UMAP with RUVr results
# mat <- ruvr_bc_norm[top_genes,] |> apply(1, scale)
# 
# pcaout <- prcomp(mat, center = F, scale. = F)
# plot(pcaout)
# pca_var<-pcaout$sdev**2
# pve <- pca_var/sum(pca_var)
# plot(pve)
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types, point.size = 4) + 
#   ggthemes::theme_par() + ggprism::scale_color_prism() + 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# um <- uwot::umap(X = pcaout$x[,1:10], n_neighbors = 5, metric = "cosine", min_dist = 0.01, spread = 5)
# prostate_meta[,c("umap_1", "umap_2")] <- um
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))
# ggplot() + 
#   geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), color = "grey", size = 4, shape = 21) + 
#   ggprism::scale_fill_prism() + 
#   ggthemes::theme_par() +
#   guides(fill = guide_legend(override.aes = list(size = 5)))

## Re-doing PCA and UMAP with limma results
mat <- limma_bc[top_genes,] |> apply(1, scale)

pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_patient_batch_corrected.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_sub_type_batch_corrected.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:10], n_neighbors = 5)
prostate_meta[,c("umap_1", "umap_2")] <- um
ggplot() +
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = patient_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_patient_batch_corrected.pdf", height = 6, width = 8)
ggplot() +
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_sub_type_batch_corrected.pdf", height = 6, width = 8)

# Hierarchical clustering
pdf(file = "prostate_heatmap_unsupervised_hclust_1500hvgs_batch_corrected.pdf", width = 12, height = 12)
rownames(mat) <- paste(prostate_meta$patient_deid, prostate_meta$tma_core_number, prostate_meta$roi, sep = "_")
ha <- HeatmapAnnotation(sub_type = prostate_meta$sub_types, col = list("sub_type" = subtype_cols))
Heatmap(matrix = mat |> t(),
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Batch\nCorrected\nExpression", 
        show_column_names = T, 
        column_names_gp = gpar(fontsize = 5),
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51))
)
dev.off()

## BLADDER (TMA2) --------------------------------------------------------------
# We will follow the same methods as we used for prostate.
bladder_meta <- meta |> filter(cancer_type == "bladder")
bladder_cts <- cts[,rownames(bladder_meta)]
bladder_norm <- norm[,rownames(bladder_meta)] 

# Colors
subtype_cols <- InSituType::colorCellTypes(names = unique(bladder_meta$sub_types), palette = "tableau20")
patient_cols <- ggprism::ggprism_data$colour_palettes$colors
names(patient_cols) <- unique(bladder_meta$patient_deid)

# HVGs, using DESeq2 to estimate dispersion# HVGs, usiprism_light2ng DESeq2 to estimate dispersion
dds <- DESeq2::DESeqDataSetFromMatrix(countData = bladder_cts[rownames(bladder_cts) != "NegProbe-WTX",], 
                                      colData = bladder_meta, design = ~1)
DESeq2::sizeFactors(dds) <- bladder_meta$q_norm_qfactors # Use the Q3 size factors
keep <- rowMeans(DESeq2::counts(dds, normalized = T)) >= 10
dds <- dds[keep,]
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- S4Vectors::mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- bladder_norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = bladder_meta$patient_deid, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_bladder_by_patient.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = bladder_meta$sub_types, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_bladder_by_sub_type.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:20], n_neighbors = 5)
bladder_meta[,c("umap_1", "umap_2")] <- um
ggplot() + 
  geom_point(data = bladder_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = patient_cols) + 
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_bladder_by_patient.pdf", height = 6, width = 8)
ggplot() + 
  geom_point(data = bladder_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_bladder_by_sub_type.pdf", height = 6, width = 8)

# There seems to be less patient influence for these samples, compared to the prostate samples. We will still use 
# batch correction for visualization.

## limma
mm <- model.matrix(~sub_types, data=bladder_meta)
limma_bc <- limma::removeBatchEffect(x = bladder_norm, batch = bladder_meta$patient_deid, design = mm)

## Finding 1500 variable genes, using non-batch-corrected data, modeling with patient, and sorting by MSE
X <- model.matrix(~patient_deid, data = bladder_meta)
hvg_resvar <- apply(bladder_norm, 1, function(y) {
  fit <- lm.fit(X, y)
  mean(fit$residuals^2)
})
top_genes <- names(sort(hvg_resvar, decreasing = TRUE))[1:1500]

## Re-doing PCA and UMAP with limma results
mat <- limma_bc[top_genes,] |> apply(1, scale)

pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = bladder_meta$patient_deid, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_bladder_by_patient_batch_corrected.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = bladder_meta$sub_types, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_bladder_by_sub_type_batch_corrected.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:10], n_neighbors = 5)
bladder_meta[,c("umap_1", "umap_2")] <- um
ggplot() +
  geom_point(data = bladder_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = patient_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_bladder_by_patient_batch_corrected.pdf", height = 6, width = 8)
ggplot() +
  geom_point(data = bladder_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_bladder_by_sub_type_batch_corrected.pdf", height = 6, width = 8)

# Hierarchical clustering
pdf(file = "bladder_heatmap_unsupervised_hclust_1500hvgs_batch_corrected.pdf", width = 12, height = 12)
rownames(mat) <- paste(bladder_meta$patient_deid, bladder_meta$tma_core_number, bladder_meta$roi, sep = "_")
ha <- HeatmapAnnotation(sub_type = bladder_meta$sub_types, col = list("sub_type" = subtype_cols))
Heatmap(matrix = mat |> t(),
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Batch\nCorrected\nExpression", 
        show_column_names = T, 
        column_names_gp = gpar(fontsize = 5),
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51))
)
dev.off()

# This looks really good. I think that using batch corrected counts for visualization from here on out will be 
# a good choice. We will always model on non-batch-corrected data though!

## Amending the metadata -------------------------------------------------------
meta <- filter(meta, sub_types != "HGPIN")
meta$sub_types_v2 <- meta$sub_types
meta$sub_types_v2[meta$cancer_type == "prostate" & meta$roi %in% c(6, 10, 13, 16, 29, 44, 47, 50)] <-  "Acinar_non_crib_LG"
meta$sub_types_v2[meta$cancer_type == "prostate" & meta$roi %in% c(2, 20, 21, 23, 26, 33, 37, 41, 54)] <-  "Acinar_non_crib_HG"
write.csv(meta, file = "meta_amended.csv")

## Redoing the plots from above, with the amended metadata ---------------------
# We will follow the same methods as we used for prostate.
prostate_meta <- meta |> filter(cancer_type == "prostate")
prostate_cts <- cts[,rownames(prostate_meta)]
prostate_norm <- norm[,rownames(prostate_meta)] 

# Colors
subtype_cols <- InSituType::colorCellTypes(names = unique(prostate_meta$sub_types_v2), palette = "tableau20")
patient_cols <- ggprism::ggprism_data$colour_palettes$colors
names(patient_cols) <- unique(prostate_meta$patient_deid)

# HVGs, using DESeq2 to estimate dispersion# HVGs, usiprism_light2ng DESeq2 to estimate dispersion
dds <- DESeq2::DESeqDataSetFromMatrix(countData = prostate_cts[rownames(prostate_cts) != "NegProbe-WTX",], 
                                      colData = prostate_meta, design = ~1)
DESeq2::sizeFactors(dds) <- prostate_meta$q_norm_qfactors # Use the Q3 size factors
keep <- rowMeans(DESeq2::counts(dds, normalized = T)) >= 10
dds <- dds[keep,]
dds <- DESeq2::estimateDispersions(dds)
DESeq2::plotDispEsts(object = dds)
disps <- S4Vectors::mcols(dds)$dispGeneEst
names(disps) <- rownames(dds)
top_genes <- sort(disps, decreasing = T)[1:1500] |> names()

# PCA
mat <- prostate_norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_patient_amended_metadata.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types_v2, point.size = 2) + 
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_sub_type_amended_metadata.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:20], n_neighbors = 5)
prostate_meta[,c("umap_1", "umap_2")] <- um
ggplot() + 
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = patient_cols) + 
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_patient_amended_metadata.pdf", height = 6, width = 8)
ggplot() + 
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types_v2), size = 2, color = "black", shape = 21) + 
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_sub_type_amended_metadata.pdf", height = 6, width = 8)

## limma
mm <- model.matrix(~sub_types_v2, data=prostate_meta)
limma_bc <- limma::removeBatchEffect(x = prostate_norm, batch = prostate_meta$patient_deid, design = mm)

## Finding 1500 variable genes, using non-batch-corrected data, modeling with patient, and sorting by MSE
X <- model.matrix(~patient_deid, data = prostate_meta)
hvg_resvar <- apply(prostate_norm, 1, function(y) {
  fit <- lm.fit(X, y)
  mean(fit$residuals^2)
})
top_genes <- names(sort(hvg_resvar, decreasing = TRUE))[1:1500]

## Re-doing PCA and UMAP with limma results
mat <- limma_bc[top_genes,] |> apply(1, scale)

pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$patient_deid, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = patient_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_patient_batch_corrected_amended_metadata.pdf", height = 6, width = 8)
ggbiplot(pcaout, var.axes = F, groups = prostate_meta$sub_types_v2, point.size = 2) +
  ggthemes::theme_par() + scale_color_manual(values = subtype_cols) +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "PCA_prostate_by_sub_type_batch_corrected_amended_metadata.pdf", height = 6, width = 8)

# UMAP
um <- uwot::umap(X = pcaout$x[,1:15], n_neighbors = 5)
prostate_meta[,c("umap_1", "umap_2")] <- um
ggplot() +
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = patient_deid), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = patient_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_patient_batch_corrected_amended_metadata.pdf", height = 6, width = 8)
ggplot() +
  geom_point(data = prostate_meta, mapping = aes(x = umap_1, y = umap_2, fill = sub_types_v2), size = 2, color = "black", shape = 21) +
  scale_fill_manual(values = subtype_cols) +
  ggthemes::theme_par() +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "UMAP_prostate_by_sub_type_batch_corrected_amended_metadata.pdf", height = 6, width = 8)

# Hierarchical clustering
pdf(file = "prostate_heatmap_unsupervised_hclust_1500hvgs_batch_corrected_amended_metadata.pdf", width = 12, height = 12)
rownames(mat) <- paste(prostate_meta$patient_deid, prostate_meta$tma_core_number, prostate_meta$roi, sep = "_")
ha <- HeatmapAnnotation(sub_type = prostate_meta$sub_types_v2, col = list("sub_type" = subtype_cols))
Heatmap(matrix = mat |> t(),
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Batch\nCorrected\nExpression", 
        show_column_names = T, 
        column_names_gp = gpar(fontsize = 5),
        show_row_names = F, 
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51))
)
dev.off()

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
#   [1] patchwork_1.3.0       magrittr_2.0.3        ggbiplot_0.6.2        ComplexHeatmap_2.22.0 ggplot2_3.5.2         dplyr_1.1.4          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4                gridExtra_2.3               rlang_1.1.5                 clue_0.3-66                
# [5] GetoptLong_1.0.5            GiottoUtils_0.2.5           matrixStats_1.5.0           compiler_4.4.2             
# [9] spatstat.geom_3.3-6         png_0.1-8                   vctrs_0.6.5                 stringr_1.5.1              
# [13] pkgconfig_2.0.3             shape_1.4.6.1               crayon_1.5.3                backports_1.5.0            
# [17] XVector_0.46.0              labeling_0.4.3              UCSC.utils_1.2.0            purrr_1.0.4                
# [21] zlibbioc_1.52.0             GenomeInfoDb_1.42.3         jsonlite_2.0.0              pak_0.9.0                  
# [25] SnowballC_0.7.1             DelayedArray_0.32.0         spatstat.utils_3.1-3        BiocParallel_1.40.2        
# [29] ggprism_1.0.6               irlba_2.3.5.1               parallel_4.4.2              cluster_2.1.6              
# [33] R6_2.6.1                    stringi_1.8.4               RColorBrewer_1.1-3          spatstat.data_3.1-6        
# [37] limma_3.62.2                reticulate_1.42.0           spatstat.univar_3.1-2       GenomicRanges_1.58.0       
# [41] Rcpp_1.0.14                 SummarizedExperiment_1.36.0 iterators_1.0.14            IRanges_2.40.1             
# [45] FNN_1.1.4.1                 Matrix_1.7-3                tidyselect_1.2.1            viridis_0.6.5              
# [49] rstudioapi_0.17.1           abind_1.4-8                 doParallel_1.0.17           codetools_0.2-20           
# [53] lattice_0.22-6              tibble_3.2.1                Biobase_2.66.0              withr_3.0.2                
# [57] askpass_1.2.1               polyclip_1.10-7             circlize_0.4.16             mclust_6.1.1               
# [61] pillar_1.10.1               lsa_0.73.3                  BiocManager_1.30.25         MatrixGenerics_1.18.1      
# [65] checkmate_2.3.2             renv_1.1.1                  foreach_1.5.2               stats4_4.4.2               
# [69] generics_0.1.3              S4Vectors_0.44.0            munsell_0.5.1               scales_1.3.0               
# [73] gtools_3.9.5                glue_1.8.0                  InSituType_2.0              tools_4.4.2                
# [77] data.table_1.16.4           RSpectra_0.16-2             locfit_1.5-9.12             umap_0.2.10.0              
# [81] colorspace_2.1-1            SingleCellExperiment_1.28.1 GenomeInfoDbData_1.2.13     cli_3.6.4                  
# [85] viridisLite_0.4.2           ggthemes_5.1.0              S4Arrays_1.6.0              uwot_0.2.3                 
# [89] gtable_0.3.6                DESeq2_1.46.0               digest_0.6.37               progressr_0.15.1           
# [93] BiocGenerics_0.52.0         SparseArray_1.6.2           rjson_0.2.23                farver_2.1.2               
# [97] lifecycle_1.0.4             httr_1.4.7                  statmod_1.5.0               GlobalOptions_0.1.2        
# [101] openssl_2.3.2              

