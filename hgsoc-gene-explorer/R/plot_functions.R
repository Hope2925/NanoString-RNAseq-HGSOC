library(ggplot2)
library(cowplot)
library(ggpubr) 
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(colorRamp2)
#library(circlize)
color_dict = list("Nano\nManso"="#FDAE6B", "Nano\nJames"="#D94801", "Nano\nBitler"="#7F2704",
        "RNAseq\nAdzib"="#3f007d", "RNAseq\nJav"="#bcbddc", "RNAseq\nEGA"="#807dba",
        "Microarray\nJimSanchez"="gold",
        "Nano Manso"="#FDAE6B", "Nano James"="#D94801", "Nano Bitler"="#7F2704",
        "RNAseq Adzib"="#3f007d", "RNAseq Jav"="#bcbddc", "RNAseq EGA"="#807dba",
        "Microarray JimSanchez"="gold",
        "PFS>12"="#008B45", "PFS<=12"="grey"
)


#########################
#### Expression #########
#########################

# Matched Gene Expression (Pre and Post) Line graphs
plot_matched_gene <- function(gene_of_interest, plot_df) {
# Plot
p <- ggplot(plot_df, aes(x = factor(Graphing_Time, levels=c("Pre", "Post", "Pre_", "Post_")), y = Value, 
                    color = factor(PFS_group, levels=c(">12mth", "<=12mth")), group = Graphing_Group)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_line(position = position_dodge(width = 0.4), aes(group = interaction(Patient, Experiment)), color = "grey50") +
  labs(title = paste("Pre vs Post for", gene_of_interest),
       y = "Expression", x = "Time", color="PFS") +
    scale_color_manual(values=c(">12mth"="springgreen4", "<=12mth"="grey")) +
    # ["RNAseq_A", "RNAseq_J", "RNAseq_EGA", "Nano_M", "Nano_B", "Nano_J"]
  theme_bw(base_size = 20) + facet_wrap(~factor(Experiment, levels=c("RNAseq Adzib", "RNAseq Jav", "RNAseq EGA", 
                                                                     "Nano Manso", "Nano Bitler", "Nano James" , 
                                                                     "Microarray JimSanchez"
                                                                     )), scale="free_y", nrow=3)

print(p)
    }

# Correlation of Experession with the PFS months
plot_exp_corr_PFS <- function(gene_of_interest, plot_df) {
  p1 <- ggplot(plot_df[plot_df$Time == "Pre",], aes(x=log10(Value+0.01), y=log10(PFS_mths+0.5))) + 
  geom_point(aes(color=Experiment)) + theme_bw(base_size=20) + 
  labs(x=paste0(gene_of_interest,": log10(Pre_Exp+0.01)"), y="log10(PFS+0.5)") + facet_wrap(~Assay, scales="free_x") + 
  geom_hline(yintercept=log10(12.01)) + scale_color_manual(values=color_dict) + 
  geom_smooth(method = "lm", se = TRUE, color="black") +  
  stat_regline_equation(size=5, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                              label.y = -0.5) + 
                facet_wrap(~Assay, scales="free_x") 

  p2 <- ggplot(plot_df[plot_df$Time == "Post",], aes(x=log10(Value+0.01), y=log10(PFS_mths+0.5))) + 
  geom_point(aes(color=Experiment)) + theme_bw(base_size=20) + 
  labs(x=paste0(gene_of_interest,": log10(Post_Exp+0.01)"), y="log10(PFS+0.5)") + 
  geom_hline(yintercept=log10(12.01)) + scale_color_manual(values=color_dict) +
  geom_smooth(method = "lm", se = TRUE, color="black") +  
  stat_regline_equation(size=5, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                              label.y = -0.5) + 
                facet_wrap(~Assay, scales="free_x") 

  print(plot_grid(p1, p2, nrow=2))
}

######################
#### L2FC ############
######################
plot_l2fc <- function(gene_of_interest, plot_l2fc_df) {
  # get the Microarray data
  plot_l2fc_micro_df = plot_l2fc_df[plot_l2fc_df$Experiment %in% c("Microarray\nJimSanchez"),]
  plot_l2fc_df = plot_l2fc_df[!plot_l2fc_df$Experiment %in% c("Microarray\nJimSanchez"),]
	# --- 1) Compute t-tests per Experiment ---
	t_test_exp <- plot_l2fc_df %>%
	group_by(Experiment) %>%
	summarise(
		t_stat  = tryCatch(t.test(Gene ~ PFS_bool)$statistic, error=function(e) NA),
		p_value = tryCatch(t.test(Gene ~ PFS_bool)$p.value,    error=function(e) NA)
	) %>%
	mutate(
		label = sprintf("p=%.2g", p_value)
	)

  t_test_micro <- plot_l2fc_micro_df %>% 
  group_by(Experiment) %>%
	summarise(
		t_stat  = tryCatch(t.test(Gene ~ PFS_bool)$statistic, error=function(e) NA),
		p_value = tryCatch(t.test(Gene ~ PFS_bool)$p.value,    error=function(e) NA)
	) %>%
	mutate(
		label = sprintf("p=%.2g", p_value)
	)

	# --- 2) Compute global t-test (all experiments pooled) ---
	t_global <- tryCatch(t.test(plot_l2fc_df$Gene ~ plot_l2fc_df$PFS_bool), error=function(e) NULL)

	global_label <- sprintf(
	 "Global p=%.2g", t_global$p.value
	)

	# --- 3) Base plot ---
	p <- ggplot(plot_l2fc_df, aes(x=Experiment, y=Gene)) +
	geom_boxplot(aes(fill=PFS_bool)) +
	geom_hline(yintercept=0) +
	labs(
		y=paste(gene_of_interest, " Log2FC"), fill="PFS Group"
	) +
	scale_fill_manual(values=c(">12mths"="springgreen4", "<=12mths"="grey")) +
	theme_bw(base_size=20)

  p_micro <- ggplot(plot_l2fc_micro_df, aes(x=Experiment, y=Gene)) +
	geom_boxplot(aes(fill=PFS_bool)) +
	geom_hline(yintercept=0) +
	labs(
		y=paste(gene_of_interest, " Log2FC"), fill="PFS Group"
	) +
	scale_fill_manual(values=c(">12mths"="springgreen4", "<=12mths"="grey")) +
	theme_bw(base_size=20)

	# --- 4) Add per-experiment t-test labels ---
	p <- p + geom_text(
	data = t_test_exp,
	aes(x = Experiment,
		y = max(plot_l2fc_df$Gene, na.rm=TRUE) * 1.05,
		label = label
	),
	size = 4.5
	)

  p_micro <- p_micro + geom_text(
	data = t_test_micro,
	aes(x = Experiment,
		y = max(plot_l2fc_micro_df$Gene, na.rm=TRUE) * 1.05,
		label = label
	),
	size = 4.5
	)

	# --- 5) Add global t-test annotation above the whole plot ---
	p <- p + annotate("text",
	x = 2.3, y = max(plot_l2fc_df$Gene, na.rm=TRUE) * 1.6,
	label = global_label, hjust = 1, vjust = 1, size = 6
	)


	#print(p)
  print(plot_grid(
  p,
  plot_grid(p_micro, NULL, rel_widths = c(0.5, 0.5)),
  ncol = 1,
  rel_heights = c(1, 1)
))
}

plot_l2fc_corr_PFS <- function(gene_of_interest, plot_l2fc_df) {
  # get the Microarray data
  plot_l2fc_micro_df = plot_l2fc_df[plot_l2fc_df$Experiment %in% c("Microarray\nJimSanchez"),]
  plot_l2fc_df = plot_l2fc_df[!plot_l2fc_df$Experiment %in% c("Microarray\nJimSanchez"),]
  # All the data
  p1 <- ggplot(plot_l2fc_df, aes(x=Gene, y=log10_PFS_mths)) + 
  geom_point(aes(color=Experiment)) + theme_bw(base_size=17) + 
  labs(x=paste(gene_of_interest, " Log2FC"), y="log10(PFS+0.5)") + #facet_wrap(~Assay, scales="free_x") + 
  geom_hline(yintercept=log10(12.01)) + geom_vline(xintercept=log2(1)) + scale_color_manual(values=color_dict) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                              label.x = min(plot_l2fc_df$Gene), label.y = -0.5, size=6)  

  p1_micro <- ggplot(plot_l2fc_micro_df, aes(x=Gene, y=log10_PFS_mths)) + 
  geom_point(aes(color=Experiment)) + theme_bw(base_size=16) + 
  labs(x=paste(gene_of_interest, " Log2FC"), y="log10(PFS+0.5)") + #facet_wrap(~Assay, scales="free_x") + 
  geom_hline(yintercept=log10(12.01)) + geom_vline(xintercept=log2(1)) + scale_color_manual(values=color_dict) + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                              label.x = min(plot_l2fc_micro_df$Gene), label.y = -0.5, size=6)  

  # Facet by Experiment
  p2 <- ggplot(plot_l2fc_df, aes(x=Gene, y=log10_PFS_mths, color=Experiment)) + 
    geom_point() + theme_bw(base_size=15) + 
    labs(x=paste(gene_of_interest, " Log2FC"), y="log10(PFS+0.5)") + #facet_wrap(~Assay, scales="free_x") + 
    geom_hline(yintercept=log10(12.01)) + geom_vline(xintercept=log2(1)) + scale_color_manual(values=color_dict) + 
    geom_smooth(method = "lm", se = TRUE) + facet_wrap(~Experiment) +
    stat_regline_equation(size=5, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")),
                                label.x = min(plot_l2fc_df$Gene), label.y = -0.5)  

  top_row <- plot_grid(
  p1, p1_micro,
  ncol = 2,
  rel_widths = c(1, 1)
)

final_plot <- plot_grid(
  top_row,
  p2,
  ncol = 1,
  rel_heights = c(1, 2)
)

print(final_plot)
}

#########################
#### Matched #########
#########################
plot_gene_scatter_helper <- function(df_use, ref_df_use, gene_use, patients, name_use, outlier_sd_use) {
	x <- log(as.numeric(df_use[gene_use, patients])+ 0.1)
    y <- log(as.numeric(ref_df_use[gene_use, patients])+0.1)
    
    data_plot <- data.frame(
      Patient = patients,
      X = x,
      Y = y
    )

	# Linear model for residuals
    fit <- lm(Y ~ X, data = data_plot)
    residuals <- resid(fit)
    sd_resid <- sd(residuals, na.rm = TRUE)
    outliers <- abs(residuals) > outlier_sd_use * sd_resid
    data_plot$Outlier <- ifelse(outliers, data_plot$Patient, NA)
	#print(data_plot[1:2,])
    
    # Correlation
    spearman_corr <- cor(data_plot$X, data_plot$Y, method = "spearman", use = "pairwise.complete.obs")
    
    p <- ggplot(data_plot, aes(x = X, y = Y)) +
      geom_point(color = "steelblue") +
      geom_smooth(method = "lm", se = FALSE, color = "red", formula = 'y ~ x') +
      geom_text_repel(aes(label = Outlier), na.rm = TRUE, color = "darkred") +
      labs(
        title = paste0("Gene: ", gene_use, " | ", name_use),
        subtitle = paste0("Spearman r = ", round(spearman_corr, 3), " (N=",length(patients) ,")"),
        x = paste0("log(", name_use, " Exp + 0.1)"),
        y = "log(NanoString Exp + 0.1)"
      ) +
      theme_bw(base_size = 14)

	  return(p)

}

plot_gene_scatter_main <- function(gene, df_list, ref_df, outlier_sd = 2) {
  # plots gene data across the different approaches with pearson correlations and outliers and fitted model
  # oultier based on sd > 2
  if (!gene %in% rownames(ref_df)) {
    stop(paste("Gene", gene, "not found in reference dataframe"))
  }

  plots <- list()
  for (name in setdiff(names(df_list), "Nano")) {
    df <- df_list[[name]]

    if (!gene %in% rownames(df)) {
      warning(paste("Gene", gene, "not found in dataframe", name))
      next
    }

    # Ensure common patients (columns)
    common_patients <- intersect(colnames(df), colnames(ref_df))
    p <- plot_gene_scatter_helper(df, ref_df, gene, common_patients, name, outlier_sd)
    plots[[name]] <- p
  }

  # Combine with patchwork if multiple plots
  if (length(plots) > 1) {
    print(Reduce(`+`, plots))
  } else {
    print(plots[[1]])
  }
}

#########################
#### SC #########
#########################

graph_gene_prepost_heatmap <- function(
    gene_mat, # should have full, pre, post with rownames=genes, and colnames=ct in ct_order
    gene,
    ct_order = c("All", "EOC_C1", "EOC_C2", "EOC_C3", "EOC_C4", "EOC_C5", "EOC_C6", "EOC_C7", "EOC_C8", "EOC_C9", 
"EOC_C10", "EOC_C11", "EOC_C12", 

"CAF-1", "CAF-2", "CAF-3", "Mesothelial", "Endothelial",

"B-cells", "Plasma-cells", "T-cells", "NK", 
"ILC", "Mast-cells", "pDC", "DC-1", "DC-2", "Macrophages")
) {

    ## COLOR FUNCTION
		col_fun <- colorRamp2(
	  c(0,  max(gene_mat, na.rm=TRUE)),
	  c("white", "red")
	  )

    ## PLOT
    ct_annot_vector = c("All", rep("EOC",12), rep("Stroma", 5), rep("Immune", 10))
    ct_anno <- HeatmapAnnotation(
      CellType = ct_annot_vector, # named vector or factor, length = nrow(mat)
      col = list(CellType = c("All"="black","EOC"="#6a6868", "Stroma"="grey20", "Immune"="grey80"))
    )
    p <- Heatmap(
        gene_mat,
        name = paste0(gene, " CPM"),
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_names_rot = 45,
		top_annotation = ct_anno,
        column_split = ct_annot_vector,
        rect_gp = gpar(col = "black", lwd = 0.5),
        row_names_side = "left", 

    )
	grid::grid.newpage()
	print(ComplexHeatmap::draw(p, newpage = TRUE))
}
