### J9/J10 Presentation Analysis ###

rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

cts <- read.csv("counts.csv", row.names = 1) |> as.matrix()
meta <- read.csv("meta.csv", row.names = 1)

all(rownames(meta) == colnames(cts)) # TRUE

library(edgeR)
library(limma)
library(dplyr)
library(magrittr)

# Just sticking with bladder for this...
bladdercts <- cts[,meta$cancer_type == "bladder"]
bladdermeta <- meta[meta$cancer_type == "bladder",]

# We will show the example of MP vs SM
bsubcts <- bladdercts[,bladdermeta$sub_types %in% c("MP", "SM")]
bsubmeta <- bladdermeta[colnames(bsubcts),]
table(bsubmeta$patient_deid, bsubmeta$sub_types)

# voom and model fitting
keep <- rowMeans(bsubcts) >= -Inf & rowMeans(bsubcts) <= 100 & rownames(bladdercts) != "NegProbe-WTX" # 11634 genes kept
cts_filt <- bsubcts[keep,]

# Size factors via TMM which is typical for sequencing data
y <- DGEList(counts = cts_filt) |> edgeR::calcNormFactors(method = "upperquartile")
mm <- model.matrix(~0+sub_types, data = bsubmeta)
vlmseq <- voomLmFit(counts = y, normalize.method = "none",
                    block = meta$Patient_number, 
                    design = mm,
                    sample.weights = T,
                    adaptive.span = T, plot = T, save.plot = T, keep.EList = T)

contrast.matrix <- makeContrasts(sub_typesSM - sub_typesMP, levels = mm)
out <- vlmseq |> contrasts.fit(contrasts = contrast.matrix) |> eBayes() |> topTable(number = Inf, coef = 1)

library(ComplexHeatmap)
toplot <- out[out$adj.P.Val < 0.05 & abs(out$logFC) > 1,] |> arrange(logFC) |> rownames()
plotidx <- arrange(bsubmeta, desc(sub_types), patient_deid) |> rownames()
mat <- vlmseq$EList$E[toplot, plotidx] |> apply(MARGIN = 1, FUN = scale) |> t()

png(filename = "8_stuff_for_presentation/MP_vs_SM_heatmap.png", width = 6, height = 10, units = "in", res = 300)
Heatmap(matrix = mat, show_row_names = F, name = "scaled\nlogCPM", cluster_rows = T, cluster_columns = T,
        # column_split = bsubmeta[plotidx,]$sub_types,
        top_annotation = HeatmapAnnotation(group = bsubmeta[plotidx,]$sub_types,
                                           patient = bsubmeta[plotidx,]$Patient_number,
                                           col = list("group" = c("SM" = "darkblue", "MP" = "limegreen")), 
                                           show_annotation_name = F, show_legend = F),
        col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), 
                                   breaks = seq(quantile(mat, 0.015), quantile(mat, 0.995), length.out = 101)))
dev.off()

library(ggplot2)
out$target <- rownames(out)
ggplot() + 
  #scattermore::geom_scattermore(data = out[(abs(out$logFC) <= 1),], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", pointsize = 2) + 
  geom_point(data = out[(abs(out$logFC) <= 1),], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", shape = ".") + 
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC > 1,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "darkblue", size = 2, shape = 16) + 
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC < -1,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "limegreen", size = 2, shape = 16) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  ggrepel::geom_text_repel(data = out[out$adj.P.Val < 0.05 & (abs(out$logFC) > 1),], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), min.segment.length = 0, box.padding = 0.25, max.overlaps = 15) +
  ggthemes::theme_par() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(filename = "8_stuff_for_presentation/MP_vs_SM_volcano.pdf", width = 10, height = 10)

mean(out$logFC > 0) # 0.6395909
mean(out$adj.P.Val < 0.05) # 0.1704487

library(msigdbr)
hallmarkgenesets <- msigdbr(db_species = "HS", collection = "H") |> 
  dplyr::select(gs_name, gene_symbol)

prlist <- arrange(out, desc(t)) |> pull(t)
names(prlist) <- arrange(out, desc(t)) |> pull(target)

library(clusterProfiler)
set.seed(2001)
hallmark <- GSEA(geneList = prlist, TERM2GENE = hallmarkgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
hallmark@result$ID %<>% factor(levels = (arrange(hallmark@result, desc(NES)) |> pull(ID)))

dat <- hallmark@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_bar(mapping = aes(y = ID, x = NES, fill = updown), stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID |> gsub(x = _, pattern = "HALLMARK_", replacement = "")),
            size = 5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID |> gsub(x = _, pattern = "HALLMARK_", replacement = "")),
            size = 5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_fill_manual(values = c("limegreen", "darkblue")) +
  ggthemes::theme_par() +
  theme(axis.text.y.left = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20), plot.title = element_text(size = 20)) +
  Seurat::NoLegend() +
  labs(title = "Hallmark Gene Sets")
ggsave(filename = "8_stuff_for_presentation/MP_vs_SM_pathways.pdf", width = 10, height = 10)

gseaplot(x = hallmark, geneSetID = "HALLMARK_G2M_CHECKPOINT", color.line = "darkblue")
gseaplot(x = hallmark, geneSetID = "HALLMARK_E2F_TARGETS", color.line = "darkblue") + 
  ggthemes::theme_par()
ggsave(filename = "8_stuff_for_presentation/MP_vs_SM_HALLMARK_E2F_TARGETS_GSEA_plot.pdf", width = 6, height = 7)
gseaplot(x = hallmark, geneSetID = "HALLMARK_ESTROGEN_RESPONSE_EARLY", color.line = "limegreen") + 
  ggthemes::theme_par()
ggsave(filename = "8_stuff_for_presentation/MP_vs_SM_HALLMARK_ESTROGEN_RESPONSE_EARLY_GSEA_plot.pdf", width = 6, height = 7)

#!#!# SHOW CONSENSUS CLASSIFICATION RESULTS #!#!#
classres <- openxlsx::read.xlsx(xlsxFile = "5_tissue_classification/consensusMIBC_classification_results.xlsx", sheet = 2, rowNames = T)
d <- classres[!classres$sub_types %in% c("Normal urothelium", "small cell"),]
d$sub_types<- factor(x = d$sub_types, levels = c("PC", "PUC", "UCIS", "CIUC", "MP", "SM"))
ggplot() + 
  geom_bar(data = d, mapping = aes(x = sub_types, fill = consenusClass), color = "black", linewidth = 0.25) + 
  ggthemes::theme_par() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_fill_brewer(type = "qual", palette = 1) +
  labs(x = "Pathological Subtype", y = "Count", fill = "Class")
ggsave(filename = "8_stuff_for_presentation/consensus_classification_results_barchart.pdf", width = 6, height = 6)

classres <- openxlsx::read.xlsx(xlsxFile = "5_tissue_classification/consensusMIBC_classification_results.xlsx", sheet = 1, rowNames = T)
consensusMIBC::plotCorrelations(xres = classres)

