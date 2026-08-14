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

# We know that we have reactome gene sets covered, but in the future, we will 
# benefit from having entrez IDs and HUGO symbols for all genes that are in 
# the WTA assay.

allgenes <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/wu-gu-prostate-bladder-data/Initial Dataset.xlsx", sheet = 2)
allgenes <- dplyr::filter(allgenes, CodeClass != "Negative") |> dplyr::select(TargetName, GeneID, HUGOSymbol)

wtagenes <- allgenes$HUGOSymbol
wtagenes <- stringr::str_split(string = wtagenes, pattern = ",")
names(wtagenes) <- allgenes$TargetName
wtagenes <- lapply(X = wtagenes, FUN = as.matrix, ncol = 1) |> lapply(as.data.frame)
wtagenes <- lapply(X = wtagenes, FUN = data.table::setnames, 1, "HUGOSymbol")
wtagenes <- dplyr::bind_rows(wtagenes, .id = "TargetName")

matched <- wtagenes[wtagenes$TargetName == wtagenes$HUGOSymbol,]
dplyr::n_distinct(matched) # 18676

# Okay, so we have a HUGO symbol for every gene. Having an entrez ID too will be
# helpful.

allgenes[allgenes$GeneID |> duplicated(),] |> nrow() # 53

allgenes$duplicatedID <-  allgenes$GeneID |> duplicated()
duped <- allgenes |> dplyr::filter(GeneID %in% allgenes[allgenes$GeneID |> duplicated(),]$GeneID) |> 
  dplyr::select(-HUGOSymbol) |> 
  dplyr::arrange(GeneID, duplicatedID)

# Let's check what we already did when we did the tissue classification.
entrezmap <- read.csv("TargetName_to_entrezid.csv", row.names = 1)
duped$entrezID <- plyr::mapvalues(x = duped$TargetName, from = entrezmap$TargetName, to = entrezmap$entrezid)

# I will run a query on NCBI
# duped$TargetName |> write.table(file = "duped.txt", row.names = F, col.names = F) 
fromncbi <- read.table("~/Downloads/entrez_query.tsv", sep = "\t", header = T)
duped$final_entrez <- plyr::mapvalues(x = duped$TargetName, from = fromncbi$Symbol, to = fromncbi$NCBI.GeneID)
duped[duped$TargetName == "H3-2",]$final_entrez <- 126961 # This one gene is weird... fixing manually.
duped$final_entrez |> duplicated() |> sum() # 0

# Ammending the entrez IDs
allgenes$entrezID <- NA
allgenes[allgenes$TargetName %in% duped$TargetName,]$entrezID <- plyr::mapvalues(x = allgenes[allgenes$TargetName %in% duped$TargetName,]$TargetName, from = duped$TargetName, to = duped$final_entrez)
allgenes[!(allgenes$TargetName %in% duped$TargetName),]$entrezID <- allgenes[!(allgenes$TargetName %in% duped$TargetName),]$GeneID
allgenes$entrezID |> dplyr::n_distinct() # 18676 

# Now, every single gene has a unique entrez ID, which is ideal.
# What if we just query every gene symbol from NCBI?
# allgenes$TargetName[1:10000] |> write.table(file = "allgenesymbols1.txt", row.names = F, col.names = F) 
# allgenes$TargetName[10001:18676] |> write.table(file = "allgenesymbols2.txt", row.names = F, col.names = F) 

query1 <- read.table("~/Downloads/allgenesymbols1_ncbi_query.tsv", sep = "\t", header = T)
query2 <- read.table("~/Downloads/allgenesymbols2_ncbi_query.tsv", sep = "\t", header = T)
allquery <- rbind(query1, query2)

allgenes$entrezID_v2 <- plyr::mapvalues(x = allgenes$TargetName, from = allquery$Symbol, to = allquery$NCBI.GeneID)
allgenes$entrezID_v2 |> dplyr::n_distinct() # 18676

# But, after checking, there are 281 genes without IDs. This is likely because 
# The HUGO symbols are for some other organism, technically, or because the HUGO
# symbol was updated. I checked this, and it seems to be the case.

# I will default to the original NanoString entrez ID.
allgenes[allgenes$entrezID_v2 %in% allgenes$TargetName,]$entrezID_v2 <- allgenes[allgenes$entrezID_v2 %in% allgenes$TargetName,]$entrezID
duplicated(allgenes$entrezID_v2) |> sum() # 0

# I personally trust my v2 version more, since I literally downloaded the updated
# data from NCBI directly.
final <- allgenes |> dplyr::rename(EntrezGeneID = entrezID_v2) |> dplyr::select(TargetName, HUGOSymbol, EntrezGeneID)

# Saving the final geomx gene names to entrex ID map:
write.table(x = final, file = "../GeoMx_TargetName_to_EnterezGeneID_map.txt", row.names = F)


