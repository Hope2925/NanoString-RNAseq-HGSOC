library(dplyr)
# ensure the gene is in all uppercase and remove any surrounding whitespace
sanitize_gene <- function(x) {
  toupper(trimws(x))
}


# Get format appropriate
reshape_for_gene <- function(mat, exp_name, gene, meta) {
  df <- as.data.frame(t(mat[gene, , drop = FALSE]))
  df$Patient <- gsub("_(Pre|Post)", "", rownames(df))
  # filter the patient
  df$Time <- ifelse(grepl("_Pre", rownames(df)), "Pre", "Post")
  df <- merge(df, meta[,c("PFS_mths", "Patient")], by.x="Patient", by.y="Patient")
  df$PFS_group <- "<=12mth"
  df[df$Patient %in% meta[meta$PFS_mths > 12,]$Patient,]$PFS_group <- ">12mth"
  df$Graphing_Group <- paste0(df$Patient, "_", df$PFS_group)
  df$Graphing_Time <- df$Time
  df[df$PFS_group == ">12mth" & df$Time == "Pre",]$Graphing_Time <- "Pre_"
  df[df$PFS_group == ">12mth" & df$Time == "Post",]$Graphing_Time <- "Post_"
  df$Experiment <- exp_name
  colnames(df)[2] <- "Value"
  return(df)
}

# Get the plot df
get_matched_plot_df <- function(gene_of_interest, exp_data, meta, nano_gene) {
  if (nano_gene) {
	mat_list <- list(exp_data$Manso, exp_data$James, exp_data$Bitler, 
                 exp_data$EGA, exp_data$Jav, exp_data$Adzib, exp_data$Jim)
	names(mat_list) <- c("Nano\nManso", "Nano\nJames", "Nano\nBitler", 
                     "RNAseq\nEGA", "RNAseq\nJav", "RNAseq\nAdzib", "Microarray\nJimSanchez")
} else {
	mat_list <- list(exp_data$EGA, exp_data$Jav, exp_data$Adzib, exp_data$Jim)
	names(mat_list) <- c("RNAseq\nEGA", "RNAseq\nJav", "RNAseq\nAdzib", "Microarray\nJimSanchez")
}
	

  plot_df <- bind_rows(lapply(names(mat_list), function(exp) {
	mat = mat_list[[exp]]
	if (!(gene_of_interest %in% rownames(mat))) return(NULL)
    reshape_for_gene(mat, exp, gene_of_interest, meta)

}))

  plot_df$Assay <- "Nano"
 plot_df[plot_df$Experiment %in% c("RNAseq\nAdzib", "RNAseq\nEGA", "RNAseq\nJav"),]$Assay <- "RNAseq"
 plot_df[plot_df$Experiment %in% c("Microarray\nJimSanchez"),]$Assay <- "Microarray"

  return(plot_df)
}
