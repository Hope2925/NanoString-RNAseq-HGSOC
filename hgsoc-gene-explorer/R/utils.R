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
get_plot_df <- function(gene_of_interest, exp_data, meta, nano_gene) {
  if (nano_gene) {
	mat_list <- list(exp_data$Manso, exp_data$James, exp_data$Bitler, 
                 exp_data$EGA, exp_data$Jav, exp_data$Adzib, exp_data$Jim)
	names(mat_list) <- c("Nano Manso", "Nano James", "Nano Bitler", 
                     "RNAseq EGA", "RNAseq Jav", "RNAseq Adzib", "Microarray JimSanchez")
} else {
	mat_list <- list(exp_data$EGA, exp_data$Jav, exp_data$Adzib, exp_data$Jim)
	names(mat_list) <- c("RNAseq EGA", "RNAseq Jav", "RNAseq Adzib", "Microarray JimSanchez")
}
	

  plot_df <- bind_rows(lapply(names(mat_list), function(exp) {
	mat = mat_list[[exp]]
	if (!(gene_of_interest %in% rownames(mat))) return(NULL)
    reshape_for_gene(mat, exp, gene_of_interest, meta)

}))

  plot_df$Assay <- "Nano"
 plot_df[plot_df$Experiment %in% c("RNAseq Adzib", "RNAseq EGA", "RNAseq Jav"),]$Assay <- "RNAseq"
 plot_df[plot_df$Experiment %in% c("Microarray JimSanchez"),]$Assay <- "Microarray"

  return(plot_df)
}

# Get the l2FC df
get_l2fc_df <- function(gene_of_interest, full_l2fc, nano_gene) {
  plot_l2fc_df = full_l2fc[,c("Experiment", "PFS_mths", "PFS_bool", gene_of_interest)]
	colnames(plot_l2fc_df) <- c("Experiment", "PFS_mths", "PFS_bool", "Gene")
	if (!nano_gene) {
	plot_l2fc_df = plot_l2fc_df[plot_l2fc_df$Experiment %in% c("RNAseq\nEGA", "RNAseq\nJav", "RNAseq\nAdzib"),]
	}
  plot_l2fc_df$log10_PFS_mths = log10(plot_l2fc_df$PFS_mths+0.5)
  plot_l2fc_df$Experiment_Name <- sub("\n", " ", plot_l2fc_df$Experiment)

  return(plot_l2fc_df)
}

# Get the SC CPM df
get_SC_df <- function(gene_of_interest, dfs, ct_order = c("All", "EOC_C1", "EOC_C2", "EOC_C3", "EOC_C4", "EOC_C5", "EOC_C6", "EOC_C7", "EOC_C8", "EOC_C9", 
"EOC_C10", "EOC_C11", "EOC_C12", 

"CAF-1", "CAF-2", "CAF-3", "Mesothelial", "Endothelial",

"B-cells", "Plasma-cells", "T-cells", "NK", 
"ILC", "Mast-cells", "pDC", "DC-1", "DC-2", "Macrophages")) {
  # CHECKS
    missing <- sapply(dfs, function(df) !(gene_of_interest %in% rownames(df)))
    if (any(missing)) {
        stop(
            "Gene missing in: ",
            paste(names(missing)[missing], collapse = ", ")
        )
    }

    # Get in proper order
    dfs <- lapply(dfs, function(df) {
        mat <- as.matrix(df)

        if (!is.null(ct_order)) {
            mat <- mat[, ct_order, drop = FALSE]
        }

        mat
    })

    ## EXTRACT GENE
    gene_mat <- do.call(
        rbind,
        lapply(dfs, function(mat) mat[gene_of_interest, , drop = FALSE])
    )

    rownames(gene_mat) <- names(dfs)
    return(gene_mat)
}