# Bladder Data
cts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("meta.csv", row.names = 1)
classres <- openxlsx::read.xlsx(xlsxFile = "5_tissue_classification/consensusMIBC_classification_results.xlsx", sheet = 1, rowNames = T)

bmeta <- meta[meta$cancer_type == "bladder",]
bmeta$MolecularClass <- plyr::mapvalues(x = bmeta$sample_id, from = rownames(classres), to = classres$consensusClass)

# Top 5 markers for each molecular class
library(limma)
library(edgeR)
bcts <- cts[rownames(cts) != "NegProbe-WTX" ,rownames(bmeta)]
y <- DGEList(counts = bcts) |> calcNormFactors(method = "upperquartile")
mm <- model.matrix(~MolecularClass, data = bmeta)
v <- voomLmFit(counts = y, 
               design = mm, 
               sample.weights = T,
               plot = T, 
               save.plot = T, 
               keep.EList = T
)

listofparams <- list()
for (MC in unique(bmeta$MolecularClass)) {
  listofparams[[MC]] <- mm[bmeta$MolecularClass == MC,] |> colMeans()
}
parameterizations <- listofparams |> dplyr::bind_rows() |> as.matrix()
rownames(parameterizations) <- names(listofparams)

outs <- list()
for (MC in unique(bmeta$MolecularClass)) {
    MCoi <- parameterizations[MC,]
    MCnoi <- (parameterizations[rownames(parameterizations) != MC,] |> colSums())/5
    contr <- (MCoi-MCnoi)
    out <- contrasts.fit(fit = v, contrasts = contr) |> eBayes() |> topTable(n = Inf)
    outs[[MC]] <- out
}

outs <- lapply(X = outs, FUN = function(df) {df$target <- rownames(df); return(df)})
markers <- dplyr::bind_rows(outs, .id = "MolecularClass") |> magrittr::set_rownames(NULL)

## - Need to have p<0.05
## - Need to have logFC>1

library(dplyr)
topmarks <- filter(markers, logFC > 0) |> 
  mutate(adj.P.Val = p.adjust(p = P.Value, method = "BH")) |> 
  filter(adj.P.Val < 0.01 & logFC > 1) |> 
  group_by(MolecularClass) |> 
  top_n(n = -10, wt = adj.P.Val) |> 
  arrange(MolecularClass, desc(logFC))
write.csv(x = topmarks, file = "MolecularClass_Top10Markers.csv", row.names = F)

# Using Claude Sonnet 4.5 and ChatGPT to make a "cheat sheet"
## - Inspired from: https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/cell-typing-cheat-sheets/
## DONE.

# Visualization
mindray_discrete <- c("#10BBB9", "#233987", "#030306", "#C7360B", "#F6BF15")
mindray_palette <- colorRampPalette(mindray_discrete)(256)
image(matrix(1:256, ncol=1), col = mindray_palette, axes = FALSE, main = "Mindray Doppler Palette")

library(ComplexHeatmap)
mat <- v$EList$E[topmarks$target, bmeta |> arrange(MolecularClass) |> pull(sample_id)] |> apply(MARGIN = 1, FUN = scale) |> t()

pdf(file = "MolecularClass_Top10Markers.pdf", width = 8, height = 10)
Heatmap(mat, 
        top_annotation = HeatmapAnnotation(Class = bmeta |> arrange(MolecularClass) |> pull(MolecularClass), 
                                           col = list("Class" = 
                                                 c("Ba/Sq" = "pink", 
                                                   "LumNS" = "orange", 
                                                   "LumP" = "gold", 
                                                   "LumU" = "limegreen", 
                                                   "NE-like" = "grey", 
                                                   "Stroma-rich" = "dodgerblue")
                                                 ),
                                           show_annotation_name = F),
        cluster_columns = F, cluster_rows = F,
        col = circlize::colorRamp2(breaks = seq(-4, 4, length.out = 101), 
                                   colors = colorRampPalette(mindray_discrete)(101)), 
        name = "Scaled\nExpression")
dev.off()

