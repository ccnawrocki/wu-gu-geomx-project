rm(list = ls())
.rs.restartR()
.libPaths()

# Ting wants GSEA done for the reactome gene sets.
# This is easy, since NanoString has pre-compiled a small database, linking 
# each of their target names to the reactome gene sets.

# For other gene sets, this can trickier, since we need to ensure that the 
# NanoString target names can accurately be mapped to the gene sets. If the 
# format is off, then genes will not be assigned correctly.

# See the learning script... I figured this problem out there.

# Here are the NanoString-curated reactome gene sets. I got this from the DSP
# control center.
reactome <- read.csv(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/reactomegenesets_for_geomx_data.csv", row.names = 1)
reactome_details <- openxlsx::read.xlsx("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/Build_78_+_NCBI_08122021_pathways_for_GeoMx_data.xlsx")

library(clusterProfiler)
library(dplyr)
library(ggplot2)

## Bladder Project GSEA --------------------------------------------------------
todo <- list("2.1 - MP_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_CIUC.csv", 
             "2.2 - MP_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_PUC.csv", 
             "2.3 - MP_vs_PC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_PC.csv", 
             "2.4 - MP_vs_SM" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_SM.csv", 
             "3.1 - PC_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_CIUC.csv", 
             "3.2 - PC_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_PUC.csv", 
             "3.3 - PC_vs_SM" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_SM.csv", 
             "4.1 - SM_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_SM_vs_CIUC.csv", 
             "4.2 - SM_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_SM_vs_PUC.csv")

for (comp in names(todo)) {
  
  # Doing the GSEA
  deres <- read.csv(file = todo[[comp]], row.names = 1)
  genelist <- dplyr::arrange(deres, desc(stat)) |> dplyr::pull(stat)
  names(genelist) <- dplyr::arrange(deres, desc(stat)) |> dplyr::pull(target)
  set.seed(9132001)
  reactres <- GSEA(geneList = genelist, TERM2GENE = reactome, 
                   eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  reactres@result$description <- plyr::mapvalues(x = reactres@result$Description, from = reactome_details$reactome_pathway_id, to = reactome_details$pathway_description, warn_missing = F)
  reactres@result$ID <- factor(x = reactres@result$description, levels = (arrange(reactres@result, desc(NES)) |> pull(description)))
  reactres@result |> select(-description) |> rename(description = ID, ID = Description) |> arrange(desc(NES)) |> write.csv(file = file.path("6_pathways/reactome_GSEA_results/bladder", paste(comp, "reactome_GSEA_results.csv", sep = "_")), row.names = F)
  
  # Plotting
  dat <- reactres@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -20, wt = p.adjust)
  dat$tosplit <- ifelse(nchar(dat$description) > 80, yes = T, no = F)
  dat$splitwhere <- ifelse(test = dat$tosplit, yes = dat$description |> stringr::str_count(pattern = " ") %/% 2, no = NA)
  dat$spaces <- ifelse(test = dat$tosplit, yes = gregexpr(" ", text = dat$description), no = NA)
  dat$splitter <- purrr::map2(.x = dat$spaces, .y = dat$splitwhere, .f = `[`) |> unlist()
  dat$description_new <- ifelse(test = dat$tosplit, 
                                yes = paste(stringr::str_sub(string = dat$description, start = 1, end = dat$splitter-1), stringr::str_sub(string = dat$description, start = dat$splitter+1, end = -1), sep = "\n"), 
                                no = dat$description)
  dat$tsize <- ifelse(test = dat$tosplit, yes = 2, no = 3)
  ggplot(data = dat) + 
    geom_bar(mapping = aes(y = ID, x = NES, fill = -log10(p.adjust)), stat = "identity") +
    geom_text(data = dat |> filter(NES < 0),
              mapping = aes(x = 0.1, y = ID, label = description_new),
              size = dat |> filter(NES < 0) |> pull(tsize), hjust = 0) +
    geom_text(data = dat |> filter(NES > 0),
              mapping = aes(x = -0.1, y = ID, label = description_new),
              size = dat |> filter(NES > 0) |> pull(tsize), hjust = 1) +
    scale_x_continuous(expand = c(0.1, 0.1)) +
    geom_vline(xintercept = 0) + 
    scale_fill_viridis_c() +
    ggthemes::theme_par() +
    ggpubr::labs_pubr() +
    theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(title = paste(gsub(x = comp, pattern = "_", replacement = " "), "Reactome Gene Sets", sep = " "))
  ggsave(filename = file.path("6_pathways/reactome_GSEA_results/bladder", paste(comp, "reactome_GSEA_results.pdf", sep = "_")), width = 12.5, height = 12)
  
  cat(comp, "done!\n", sep = " ")
}


## Prostate Project GSEA -------------------------------------------------------

# Same as above for all requests, except for 2.3:

todo <- list("1.1 - Acinar_non_crib_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_revised_results/07/07 - DESeq2_Acinar_non_crib_vs_Acinar_crib.csv", 
             "1.2 - Acinar_non_crib_HG_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_revised_results/07/07 - DESeq2_Acinar_non_crib_HG_vs_Acinar_crib.csv",
             "1.3 - Acinar_non_crib_LG_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_revised_results/07/07 - DESeq2_Acinar_non_crib_LG_vs_Acinar_crib.csv",
             "2.1.1 - Acinar_IDC-P_vs_all_acinar_carcinoma" = "4_DE_analysis/DE_analysis_revised_results/14/14.1 - DESeq2_all_acinar_IDC-P_vs_all_acinar_carcinoma.csv",
             "2.1.2 - Acinar_IDC-P_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_revised_results/14/14.4 - DESeq2_all_acinar_IDC-P_vs_Acinar_crib.csv",
             "2.1.3 - Acinar_IDC-P_vs_Acinar_non_crib_HG" = "4_DE_analysis/DE_analysis_revised_results/14/14.3 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_HG.csv",
             "2.1.4 - Acinar_IDC-P_vs_Acinar_non_crib_LG" = "4_DE_analysis/DE_analysis_revised_results/14/14.2 - DESeq2_all_acinar_IDC-P_vs_Acinar_non_crib_LG.csv",
             "2.2.1 - Paired_Acinar_IDC-P_crib_vs_all_acinar_carcinoma" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/13 - DESeq2_all_acinar_carcinoma_vs_Acinar_IDC-P_crib_core_paired.csv",
             "2.2.2 - Paired_Acinar_IDC-P_crib_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/13 - DESeq2_Acinar_crib_vs_Acinar_IDC-P_crib_core_paired.csv",
             "2.2.3 - Paired_Acinar_IDC-P_crib_vs_Acinar_non_crib_HG" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/13 - DESeq2_Acinar_non_crib_HG_vs_Acinar_IDC-P_crib_core_paired.csv",
             "2.2.4 - Paired_Acinar_IDC-P_crib_vs_Acinar_non_crib_LG" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/13 - DESeq2_Acinar_non_crib_LG_vs_Acinar_IDC-P_crib_core_paired.csv",
             "2.4.1 - Acinar_AIP_vs_all_acinar_carcinoma" = "4_DE_analysis/DE_analysis_revised_results/15/15.1 - DESeq2_AIP_vs_all_acinar_carcinoma.csv",
             "2.4.2 - Acinar_AIP_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_revised_results/15/15.4 - DESeq2_AIP_vs_Acinar_crib.csv",
             "2.4.3 - Acinar_AIP_vs_Acinar_non_crib_HG" = "4_DE_analysis/DE_analysis_revised_results/15/15.3 - DESeq2_AIP_vs_Acinar_non_crib_HG.csv",
             "2.4.4 - Acinar_AIP_vs_Acinar_non_crib_LG" = "4_DE_analysis/DE_analysis_revised_results/15/15.2 - DESeq2_AIP_vs_Acinar_non_crib_LG.csv",
             "2.5 - IDC-P_vs_AIP" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/16 - DESeq2_all_acinar_IDC-P_vs_AIP.csv",
             "3.1.1 - Ductal_vs_Acinar_crib" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/03 - DESeq2_Ductal_vs_Acinar_crib.csv",
             "3.1.2 - Ductal_vs_Acinar_non_crib" = "4_DE_analysis/DE_analysis_revised_results/08/08 - DESeq2_Ductal_vs_Acinar_non_crib.csv",
             "3.1.3 - Ductal_vs_Acinar_non_crib_HG" = "4_DE_analysis/DE_analysis_revised_results/08/08 - DESeq2_Ductal_vs_Acinar_non_crib_HG.csv",
             "3.1.4 - Ductal_vs_Acinar_non_crib_LG" = "4_DE_analysis/DE_analysis_revised_results/08/08 - DESeq2_Ductal_vs_Acinar_non_crib_LG.csv",
             "3.2.1 - Paired_Ductal_vs_Acinar" = "4_DE_analysis/DE_analysis_revised_results/10/10.4 - DESeq2_Ductal_vs_Acinar_core_paired.csv",
             "3.2.2 - Paired_Ductal_vs_Acinar_non_crib" = "4_DE_analysis/DE_analysis_revised_results/10/10.1 - DESeq2_Ductal_vs_Acinar_non_crib_core_paired.csv",
             "3.2.3 - Paired_Ductal_vs_Acinar_non_crib_HG" = "4_DE_analysis/DE_analysis_revised_results/10/10.3 - DESeq2_Ductal_vs_Acinar_HG_core_paired.csv",
             "3.2.4 - Paired_Ductal_vs_Acinar_non_crib_LG" = "4_DE_analysis/DE_analysis_revised_results/10/10.2 - DESeq2_Ductal_vs_Acinar_non_crib_LG_core_paired.csv",
             "4.1 - Ductal_vs_Ductal_IDC-P" = "4_DE_analysis/DE_analysis_revised_results/12/12 - DESeq2_Ductal_vs_Ductal_IDC-P.csv",
             "4.2 - Paired_Ductal_vs_Ductal_IDC-P" = "4_DE_analysis/DE_analysis_revised_results/13/13 - DESeq2_Ductal_vs_Ductal_IDC-P_core_paired.csv")

# One little issue we ran into is that a few gene sets have a couple versions, so
# there are duplicated descriptions.

# We will take the versions with the highest coverage.
reactome_details$coverage_prop <- (gsub(pattern = "\\%", replacement = "", x = reactome_details$pathway_coverage_pct) |> as.numeric())/100
reactome_details <- reactome_details |> group_by(pathway_description) |> mutate(highest_coverage = max(coverage_prop))
reactome_details <- reactome_details[reactome_details$coverage_prop == reactome_details$highest_coverage, ]

# For some reason, there are a few duplicated rows here. We will remove them.
reactome_details <- reactome_details[!duplicated(reactome_details$pathway_description),]
reactome <- reactome[reactome$gs_name %in% reactome_details$reactome_pathway_id,]

for (comp in names(todo)) {
  
  # Doing the GSEA
  deres <- read.csv(file = todo[[comp]], row.names = 1)
  genelist <- dplyr::arrange(deres, desc(stat)) |> dplyr::pull(stat)
  names(genelist) <- dplyr::arrange(deres, desc(stat)) |> dplyr::pull(target)
  set.seed(9132001)
  reactres <- GSEA(geneList = genelist, TERM2GENE = reactome, 
                   eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  reactres@result$description <- plyr::mapvalues(x = reactres@result$Description, from = reactome_details$reactome_pathway_id, to = reactome_details$pathway_description, warn_missing = F)
  reactres@result$ID <- factor(x = reactres@result$description, levels = (arrange(reactres@result, desc(NES)) |> pull(description)))
  reactres@result |> select(-description) |> rename(description = ID, ID = Description) |> arrange(desc(NES)) |> write.csv(file = file.path("6_pathways/reactome_GSEA_results/prostate", paste(comp, "reactome_GSEA_results.csv", sep = "_")), row.names = F)
  
  # Plotting
  dat <- reactres@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -20, wt = p.adjust)
  dat$tosplit <- ifelse(nchar(dat$description) > 50, yes = T, no = F)
  dat$splitwhere <- ifelse(test = dat$tosplit, yes = dat$description |> stringr::str_count(pattern = " ") %/% 2, no = NA)
  dat$spaces <- ifelse(test = dat$tosplit, yes = gregexpr(" ", text = dat$description), no = NA)
  dat$splitter <- purrr::map2(.x = dat$spaces, .y = dat$splitwhere, .f = `[`) |> unlist()
  dat$description_new <- ifelse(test = dat$tosplit, 
                                yes = paste(stringr::str_sub(string = dat$description, start = 1, end = dat$splitter-1), stringr::str_sub(string = dat$description, start = dat$splitter+1, end = -1), sep = "\n"), 
                                no = dat$description)
  dat$tsize <- ifelse(test = dat$tosplit, yes = 2, no = 3)
  ggplot(data = dat) + 
    geom_bar(mapping = aes(y = ID, x = NES, fill = -log10(p.adjust)), stat = "identity") +
    geom_text(data = dat |> filter(NES < 0),
              mapping = aes(x = 0.1, y = ID, label = description_new),
              size = dat |> filter(NES < 0) |> pull(tsize), hjust = 0) +
    geom_text(data = dat |> filter(NES > 0),
              mapping = aes(x = -0.1, y = ID, label = description_new),
              size = dat |> filter(NES > 0) |> pull(tsize), hjust = 1) +
    scale_x_continuous(expand = c(0.05, 0.05)) +
    geom_vline(xintercept = 0) + 
    scale_fill_viridis_c() +
    ggthemes::theme_par() +
    ggpubr::labs_pubr() +
    theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(title = paste(gsub(x = comp, pattern = "_", replacement = " "), "Reactome Gene Sets", sep = " "))
  ggsave(filename = file.path("6_pathways/reactome_GSEA_results/prostate", paste(comp, "reactome_GSEA_results.pdf", sep = "_")), width = 12.5, height = 12)
  
  cat(comp, "done!\n", sep = " ")
}

# 2.3: ssGSEA, due to n=1, + a revision to the DE
# Regarding statistics of Q11, you mentioned that “I can revise the heat map, but
# we should NOT do the testing or report p-values for each gene. Having n=1 on 
# one side of the comparison does not permit testing. Instead, we can show the 
# genes with the largest effect sizes, then maybe do a pathway analysis, like 
# ssGSEA.” “For point 11, I am showing genes with |log2FC| > 2. We cannot do 
# inference for this, as one side of the comparison is n=1. We can take the list 
# of genes, rank by log2FC, then do ssGSEA.”

# I agree with you. Please do ssGSEA as you suggested.

library(ssGSEA2)

norm <- read.csv("log2plus1_q3norm.csv", row.names = 1)
meta <- read.csv(file = "meta_amended.csv", row.names = 1)
norm <- norm[rownames(norm) != "NegProbe-WTX", rownames(meta)]

prostate <- meta[(meta$cancer_type == "prostate"),] |> rownames()
norm <- norm[,prostate]
meta <- meta[prostate,]

meta$sub_types_v5 <- case_when(meta$roi %in% c(48, 49) ~ "Precursor", 
                               meta$sub_types_v2 %in% c("Acinar IDC-P_crib", "AIP") ~ "Intraductal_spread", 
                               T ~ meta$sub_types_v2)

small_meta <- meta |> filter(sub_types_v5 %in% c("Precursor", "Intraductal_spread"))
small_meta$group <- small_meta$sub_types_v5

# gct <- cbind(Name = rownames(norm[,rownames(small_meta)]), Description = "na", norm[,rownames(small_meta)])
# con <- file("6_pathways/reactome_ssGSEA_results/data_for_2.3.gct", "wt")
# writeLines("#1.2", con)
# writeLines(paste(nrow(norm[,rownames(small_meta)]), ncol(norm[,rownames(small_meta)]), sep="\t"), con)
# write.table(gct, con, sep="\t", quote=FALSE, row.names=FALSE)
# close(con)
# 
# reactome_details_gmt <- select(reactome_details, pathway_description, targets)
# reactome_details_gmt$source <- "reactome"
# reactome_details_gmt$targets <- gsub(pattern = ",", replacement = "\t", x = reactome_details_gmt$targets)
# reactome_details_gmt$pathway_description <- gsub(pattern = " ", replacement = "_", x = reactome_details_gmt$pathway_description)
# reactome_details_gmt$gmt <- paste(reactome_details_gmt$pathway_description, reactome_details_gmt$source, reactome_details_gmt$targets, sep = "\t")
# write.table(reactome_details_gmt[,"gmt"], file = "6_pathways/reactome_ssGSEA_results/reactome.gmt", row.names = F, col.names = F)

# # This took ~4 hours, I used the same parameters in the GitHub vignette
# ssGSEA_res <- ssGSEA2::run_ssGSEA2(
#   input.ds = "6_pathways/reactome_ssGSEA_results/data_for_2.3.gct", 
#   gene.set.databases = "6_pathways/reactome_ssGSEA_results/reactome.gmt", 
#   output.prefix = "ssGSEA_results_for_2.3", 
#   output.directory = "6_pathways/reactome_ssGSEA_results", 
#   weight = 0.75, 
#   sample.norm.type = "none", 
#   output.score.type = "NES", 
#   correl.type = "rank", 
#   statistic = "area.under.RES",
#   nperm = 1000, 
#   min.overlap = 5, 
#   extended.output = TRUE, 
#   global.fdr = FALSE
#   )
# saveRDS(object = ssGSEA_res, file = "6_pathways/reactome_ssGSEA_results/results_for_2.3.RDS")

# Once I had the results, I deleted the files I made above. There is no point in
# keeping them... they are just duplicated data and can be remade easily.

ssGSEA_res <- readRDS("6_pathways/reactome_ssGSEA_results/results_for_2.3.RDS")
names(ssGSEA_res) <- gsub(pattern = "\"", replacement = "", x = names(ssGSEA_res))
names(ssGSEA_res) <- gsub(pattern = "_", replacement = " ", x = names(ssGSEA_res))
resmat <- lapply(ssGSEA_res, FUN = sapply, `[[`, "ES") |> 
  dplyr::bind_rows(.id = "pathway")

byvar <- (apply(X = resmat, MARGIN = 2, FUN = var) |> sort(decreasing = T))

mm <- model.matrix(~0+group, data = small_meta)
groupmeans <- (t(as.matrix(resmat)) %*% mm) |> sweep(x = _, MARGIN = 2, STATS = colSums(mm), FUN = "/")
colnames(groupmeans) <- gsub(pattern = "group", replacement = "", x = colnames(groupmeans))
groupmeans <- as.data.frame(groupmeans)
bymeandiff <- ((groupmeans[,"Precursor", drop = F] - groupmeans[,"Intraductal_spread", drop = F])) |> data.table::setnames("meandiff")
bymeandiff$pathway <- rownames(bymeandiff)
topmeandiffs <- bymeandiff |> 
  mutate(updown = sign(meandiff)) |>
  group_by(updown) |>
  top_n(n = 50, wt = abs(meandiff)) |> 
  pull(pathway)
bymeandiff$absmeandiff <- abs(bymeandiff$meandiff)
topabsmeandiffs <- bymeandiff |> arrange(desc(absmeandiff)) |> pull(pathway)

summary_table <- bymeandiff |> select(pathway, meandiff, absmeandiff)
summary_table$variance <- plyr::mapvalues(x = summary_table$pathway, from = names(byvar), to = byvar) |> as.numeric()
summary_table <- dplyr::mutate(summary_table, updown = sign(meandiff)) |> 
  dplyr::arrange(desc(updown), desc(abs(meandiff))) |> as.data.frame() |> select(-updown)
write.csv(x = summary_table, file = "2.3 - Precursor_vs_Intraductal_spread_reactome_ssGSEA_summary_table.csv", row.names = F)

resmat <- apply(X = resmat, MARGIN = 2, FUN = scale) |> t()
colnames(resmat) <- sprintf("%03d", small_meta$roi)
write.csv(x = resmat, file = "2.3 - Precursor_vs_Intraductal_spread_reactome_ssGSEA_row-scaled_enrichment_scores.csv")

library(ComplexHeatmap)
mat <- resmat[topmeandiffs,]
ha <- HeatmapAnnotation(group = small_meta$group, 
                        sub_type = small_meta$sub_types_v2,
                        col = list("group"=c("Precursor"="dodgerblue", "Intraductal_spread"="gold3"), 
                                   "sub_type"=c("Acinar IDC-P_crib"="blue", "AIP"="purple")))
pdf("2.3 - Precursor_vs_Intraductal_spread_reactome_ssGSEA_top_50_pathways_by_mean_difference_in_each_direction.pdf", width = 8, height = 16)
Heatmap(matrix = mat, 
        column_split = small_meta$group, 
        top_annotation = ha, 
        heatmap_legend_param = list(legend_direction = "horizontal"),
        row_names_gp = gpar(fontsize = 4), 
        column_names_gp = gpar(fontsize = 7),
        name = "row-scaled\nenrichment score",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
) |> 
  draw(heatmap_legend_side = "bottom", annotation_legend_side = "top")
dev.off()

mat <- resmat[topabsmeandiffs[1:100],]
pdf("2.3 - Precursor_vs_Intraductal_spread_reactome_ssGSEA_top_100_pathways_by_absolute_mean_difference.pdf", width = 8, height = 16)
Heatmap(matrix = mat, 
        column_split = small_meta$group, 
        top_annotation = ha, 
        heatmap_legend_param = list(legend_direction = "horizontal"),
        row_names_gp = gpar(fontsize = 4), 
        column_names_gp = gpar(fontsize = 7),
        name = "row-scaled\nenrichment score",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51)), 
        width = ncol(mat)*unit(3, "mm"),
        height = nrow(mat)*unit(2, "mm"), 
        column_title_gp = gpar(col = NA)
) |> 
  draw(heatmap_legend_side = "bottom", annotation_legend_side = "top")
dev.off()

pdf("2.3 - Precursor_vs_Intraductal_spread_reactome_ssGSEA_top_500_pathways_by_variance.pdf", width = 8, height = 16)
mat <- resmat[names(byvar)[1:500],]
Heatmap(matrix = mat, 
        top_annotation = ha,
        show_row_names = F,
        column_names_gp = gpar(fontsize = 7),
        name = "row-scaled\nenrichment score",
        col = circlize::colorRamp2(breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out=51),
                                   colors = viridis::viridis(51))
)
dev.off()

