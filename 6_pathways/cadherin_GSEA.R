rm(list = ls())
.rs.restartR()
.libPaths()

# Could you please do pathway analysis for specific gene sets relating to 
# cadherins for following comparisons? Thank you!
#   PC_vs_CIUC
#   PC_vs_PUC

# To do this, I will simply load every msigdb pathway. Then, I will subset down
# to only the pathways with "cadherin" in their description. This will be a 
# decent starting point.

library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)

genesets <- msigdbr(db_species = "HS")
cadherinsets <- genesets[grepl(pattern = "(cadherin)|(Cadherin)|(CADHERIN)", x = genesets$gs_description) | 
                           grepl(pattern = "(cadherin)|(Cadherin)|(CADHERIN)", x = genesets$gs_name),
                         ]
cadherinsets$gs_name |> unique()
# [1] "GOBP_ADHERENS_JUNCTION_ASSEMBLY"                                     "GOBP_ADHERENS_JUNCTION_MAINTENANCE"                                 
# [3] "GOBP_ADHERENS_JUNCTION_ORGANIZATION"                                 "GOBP_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN"                       
# [5] "GOBP_NEGATIVE_REGULATION_OF_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN" "GOBP_POSITIVE_REGULATION_OF_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN"
# [7] "GOBP_REGULATION_OF_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN"          "GOCC_ADHERENS_JUNCTION"                                             
# [9] "GOCC_CATENIN_COMPLEX"                                                "GOCC_DESMOSOME"                                                     
# [11] "GOCC_FASCIA_ADHERENS"                                                "GOMF_CADHERIN_BINDING"                                              
# [13] "GOMF_CADHERIN_BINDING_INVOLVED_IN_CELL_CELL_ADHESION"                "PID_ECADHERIN_KERATINOCYTE_PATHWAY"                                 
# [15] "PID_ECADHERIN_NASCENT_AJ_PATHWAY"                                    "PID_ECADHERIN_STABILIZATION_PATHWAY"                                
# [17] "PID_NCADHERIN_PATHWAY" 

# These are the gene sets that we will test. Now, we will get our entrez IDs, 
# which will reduce gene dropout when we run GSEA.
gs_info <- select(cadherinsets, 8:14, 17:18) |> filter(!duplicated(x = gs_name))
write.csv(gs_info, "selected_cadherin_gene_sets.csv")
cadherinsets <- select(cadherinsets, gs_id, ncbi_gene)
symboltoentrez <- read.csv("TargetName_to_entrezid.csv", row.names = 1)

##### PC vs CIUC ---------------------------------------------------------------
pc_vs_ciuc <- read.csv("4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_CIUC.csv", row.names = 1)
pc_vs_ciuc$entrezid <- plyr::mapvalues(x = pc_vs_ciuc$target, from = symboltoentrez$TargetName, to = symboltoentrez$entrezid)

genelist <- dplyr::arrange(pc_vs_ciuc, desc(stat)) |> dplyr::pull(stat)
names(genelist) <- dplyr::arrange(pc_vs_ciuc, desc(stat)) |> dplyr::pull(entrezid)
set.seed(9132001)
pc_vs_ciuc_res <- GSEA(geneList = genelist, TERM2GENE = cadherinsets, 
                       eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
pc_vs_ciuc_res@result$description <- plyr::mapvalues(x = pc_vs_ciuc_res@result$Description, from = gs_info$gs_id, to = gs_info$gs_name, warn_missing = F)
pc_vs_ciuc_res@result$ID <- factor(x = pc_vs_ciuc_res@result$description, levels = (arrange(pc_vs_ciuc_res@result, desc(NES)) |> pull(description)))

# Plotting
dat <- pc_vs_ciuc_res@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -20, wt = p.adjust)
ggplot(data = dat) + 
  geom_bar(mapping = aes(y = ID, x = NES, fill = -log10(p.adjust)), stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 3, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 3, hjust = 1) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  geom_vline(xintercept = 0) + 
  scale_fill_viridis_c() +
  ggthemes::theme_par() +
  ggpubr::labs_pubr() +
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(title = "PC vs. CIUC: Cadherin-related Gene Sets")
ggsave(filename = "PC_vs_CIUC_Cadherin-related_gene_sets.pdf", width = 8, height = 6)

gseaplot(x = pc_vs_ciuc_res, geneSetID = "M17581", title = "GOCC_DESMOSOME")
ggsave(filename = "PC_vs_CIUC_GOCC_DESMOSOME_GSEA_plot.pdf", height = 8, width = 5)
gseaplot(x = pc_vs_ciuc_res, geneSetID = "M19119", title = "GOMF_CADHERIN_BINDING")
ggsave(filename = "PC_vs_CIUC_GOMF_CADHERIN_BINDING_GSEA_plot.pdf", height = 8, width = 5)

##### PC vs PUC ----------------------------------------------------------------
pc_vs_puc <- read.csv("4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_PUC.csv", row.names = 1)
pc_vs_puc$entrezid <- plyr::mapvalues(x = pc_vs_puc$target, from = symboltoentrez$TargetName, to = symboltoentrez$entrezid)

genelist <- dplyr::arrange(pc_vs_puc, desc(stat)) |> dplyr::pull(stat)
names(genelist) <- dplyr::arrange(pc_vs_puc, desc(stat)) |> dplyr::pull(entrezid)
set.seed(9132001)
pc_vs_puc_res <- GSEA(geneList = genelist, TERM2GENE = cadherinsets, 
                       eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
pc_vs_puc_res@result$description <- plyr::mapvalues(x = pc_vs_puc_res@result$Description, from = gs_info$gs_id, to = gs_info$gs_name, warn_missing = F)
pc_vs_puc_res@result$ID <- factor(x = pc_vs_puc_res@result$description, levels = (arrange(pc_vs_puc_res@result, desc(NES)) |> pull(description)))

# Plotting
gseaplot(x = pc_vs_puc_res, geneSetID = "M19119", title = "GOMF_CADHERIN_BINDING")
ggsave(filename = "PC_vs_PUC_GOMF_CADHERIN_BINDING_GSEA_plot.pdf", height = 8, width = 5)

