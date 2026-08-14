rm(list = ls())

reslist <- list.files("4_DE_analysis/DE_analysis_second_revision_and_bladder_results", pattern = "^bladder(.+)csv$")
names(reslist) <- gsub("(bladder - DESeq2_)|(\\.csv)", "", reslist)

cdh1res <- list()
for (comp in names(reslist)) {
  tmp <- data.table::fread(file.path("4_DE_analysis/DE_analysis_second_revision_and_bladder_results", reslist[[comp]]))
  cdh1res[[comp]] <- tmp[tmp$target == "CDH1", -1]
}

cdh1res <- dplyr::bind_rows(cdh1res, .id = "comparison")
write.csv(cdh1res, file = "4_DE_analysis/bladder_results_and_CDH1.csv", row.names = F)
