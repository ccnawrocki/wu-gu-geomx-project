bladder_res_list <- 
  list("MP_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_CIUC.csv", 
    "PC_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_CIUC.csv", 
    "SM_vs_CIUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_SM_vs_CIUC.csv", 
    "MP_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_PUC.csv", 
    "SM_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_SM_vs_PUC.csv", 
    "PC_vs_PUC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_PUC.csv",
    "MP_vs_PC" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_PC.csv", 
    "MP_vs_SM" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_MP_vs_SM.csv", 
    "PC_vs_SM" = "4_DE_analysis/DE_analysis_second_revision_and_bladder_results/bladder - DESeq2_PC_vs_SM.csv"
    )
bladder_res_list <- lapply(X = bladder_res_list, FUN = read.csv, row.names = 1)
bladder_res_list <- lapply(X = bladder_res_list, FUN = function(x) {x[abs(x$log2FoldChange) > 1 & x$padj < 0.05,]})
bladder_res_list <- lapply(X = bladder_res_list, FUN = function(x) {x |> dplyr::arrange(desc(log2FoldChange)) |> as.data.frame()})

WB <- openxlsx::buildWorkbook(x = bladder_res_list, asTable = T)
openxlsx::saveWorkbook(wb = WB, file = "bladder_de_results_summarized.xlsx", overwrite = T)

