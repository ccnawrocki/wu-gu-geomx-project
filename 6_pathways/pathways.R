library(msigdbr)
hallmarkgenesets <- msigdbr(species = "Homo sapiens", category = "H") |> 
  dplyr::select(gs_name, gene_symbol)

head(hallmarkgenesets)
#                 gs_name gene_symbol
# 1 HALLMARK_ADIPOGENESIS       ABCA1
# 2 HALLMARK_ADIPOGENESIS       ABCB8
# 3 HALLMARK_ADIPOGENESIS       ACAA2
# 4 HALLMARK_ADIPOGENESIS       ACADL
# 5 HALLMARK_ADIPOGENESIS       ACADM
# 6 HALLMARK_ADIPOGENESIS       ACADS

reactome <- openxlsx::read.xlsx("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/Build_78_+_NCBI_08122021_pathways_for_GeoMx_data.xlsx")
genelist <- reactome$targets |> stringr::str_split(pattern = ",")
names(genelist) <- reactome$reactome_pathway_id
reactomegenesets <- purrr::map(.x = genelist, .f = matrix, ncol = 1, dimnames = list(NULL, "gene_symbol")) |> 
  purrr::map(.f = as.data.frame) |> 
  dplyr::bind_rows(.id = "gs_name")
# write.csv(x = reactomegenesets, file = "reactomegenesets_for_geomx_data.csv")

library(clusterProfiler)
test <- read.csv(file = "4_DE_analysis/DESeq2_Acinar_non_crib_vs_Acinar_crib.csv", row.names = 1)
testlist <- dplyr::arrange(test, desc(stat)) |> dplyr::pull(stat)
names(testlist) <- dplyr::arrange(test, desc(stat)) |> dplyr::pull(target)
reactres <- GSEA(geneList = testlist, TERM2GENE = reactomegenesets, 
                 eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
reactres@result$description <- plyr::mapvalues(x = reactres@result$Description, from = reactome$reactome_pathway_id, to = reactome$pathway_description)

library(dplyr)
library(ggplot2)
reactres@result$ID <- factor(x = reactres@result$description, levels = (arrange(reactres@result, desc(NES)) |> pull(description)))
dat <- reactres@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -20, wt = p.adjust)
ggplot(data = dat) + 
  geom_bar(mapping = aes(y = ID, x = NES, fill = -log10(p.adjust)), stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_fill_viridis_c() +
  ggthemes::theme_par() +
  ggpubr::labs_pubr() +
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "Reactome Gene Sets") 


