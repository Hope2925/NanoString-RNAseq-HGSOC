# Compare 

#######################
## HELPER FUNCTIONS ###
#######################
# Get the top X neighbors for a given gene from the TOM
get_top_neighbors <- function(tom_mat, gene, top_k = 50, min_weight = 0.05) {
  if (!gene %in% rownames(tom_mat)) return(character(0))
  v <- tom_mat[gene, ]
  v <- v[!is.na(v)]
  v <- v[names(v) != gene]
  
  if (!is.null(min_weight)) {
    keep <- names(v)[v >= min_weight]
    return(keep)
  } else {
    ord <- order(v, decreasing = TRUE)
    return(names(v)[ord][seq_len(min(top_k, length(ord)))])
  }
}

# Calculate the jaccard index
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

# Get the correlation between two vectors while accounting for NA problems
safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

# Get the number of genes that a gene is connected to
row_connectivity <- function(tom_mat, gene) {
  if (!gene %in% rownames(tom_mat)) return(NA_real_)
  v <- tom_mat[gene, ]
  v <- v[is.finite(v)]
  v <- v[names(v) != gene]
  sum(v)
}

# Use this to ensure identical gene universe between TOMs
restrict_to_shared_genes <- function(tom1, tom2, gene) {
  g <- intersect(colnames(tom1), colnames(tom2))
  cat("\nRestricting to ", length(g), " genes for ", gene)
  if (gene %in% rownames(tom1)) {
    tom1 <- tom1[gene, g, drop = FALSE]
  }
  if (gene %in% rownames(tom2)) {
    tom2 <- tom2[gene, g, drop = FALSE]
  }
  list(tom1 = tom1, tom2 = tom2)
}

# Get the proper matrix with the gene if multiple 
find_gene_block <- function(tom_list, gene) {
  for (tom in tom_list) {
    if (gene %in% rownames(tom)) {
      return(tom)
    }
  }
}

#######################
## PREP INPUT ###
#######################
genes_of_interest = c("IDO1", "CD40", "JAK2", "CTSS", "ALDOC", "TMUB2", "ID4", "ENO1", "EGR1", "OAS3", 
              "RRM2", "GBP4", "AMOTL2", "ESYT3", "HCN3", "DDX11", "GPR173", "CLK2", "CDKL2", "SHROOM1", "CD274")
length(genes_of_interest)

genen_dir <- "/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/gene_networks"
TOM_pre1 <- readRDS(file.path(genen_dir, "TOM_filt1_PRE.tsv"))
TOM_pre2 <- readRDS(file.path(genen_dir, "TOM_filt2_PRE.tsv"))
dim(TOM_pre1)
dim(TOM_pre2)

TOM_post1 <- readRDS(file.path(genen_dir, "TOM_filt1_POST.tsv"))
TOM_post2 <- readRDS(file.path(genen_dir, "TOM_filt2_POST.tsv"))
dim(TOM_post1)
dim(TOM_post2)

TOM_bulk <- list(
  Pre = list(block1 = TOM_pre1, block2 = TOM_pre2),
  Post = list(block1 = TOM_post1, block2 = TOM_post2)
)



TOM_sc <- list(
  Macrophages = list(Pre = "", Post = ""),
  T_cells     = list(Pre = "", Post = ""),
  EOC     = list(Pre = "", Post = ""), 
  CAF     = list(Pre = "", Post = ""), 
  Mesothelial = list(Pre = "", Post = ""), 
  Innate_Lymphoid = list(Pre = "", Post = ""),  # "NK", "ILC"
  B_cell = list(Pre = "", Post = ""), 
  Myeloid_APC = list(Pre = "", Post = "") # All DCs
)

TOM_sc$Macrophages$Pre <- readRDS(file.path(genen_dir, "Macrophage_Pre_Filt_TOM.rds"))
TOM_sc$Macrophages$Post <- readRDS(file.path(genen_dir, "Macrophage_Post_Filt_TOM.rds"))
TOM_sc$T_cells$Pre <- readRDS(file.path(genen_dir, "T-cells_Pre_Filt_TOM.rds"))
TOM_sc$T_cells$Post <- readRDS(file.path(genen_dir, "T-cells_Post_Filt_TOM.rds"))
TOM_sc$EOC$Pre <- readRDS(file.path(genen_dir, "EOC_Pre_Filt_TOM.rds"))
TOM_sc$EOC$Post <- readRDS(file.path(genen_dir, "EOC_Post_Filt_TOM.rds"))
TOM_sc$CAF$Pre <- readRDS(file.path(genen_dir, "CAF_Pre_Filt_TOM.rds"))
TOM_sc$CAF$Post <- readRDS(file.path(genen_dir, "CAF_Post_Filt_TOM.rds"))
TOM_sc$Mesothelial$Pre <- readRDS(file.path(genen_dir, "Mesothelial_Pre_Filt_TOM.rds"))
TOM_sc$Mesothelial$Post <- readRDS(file.path(genen_dir, "Mesothelial_Post_Filt_TOM.rds"))
TOM_sc$Innate_Lymphoid$Pre <- readRDS(file.path(genen_dir, "Innate_Lymphoid_Pre_Filt_TOM.rds"))
TOM_sc$Innate_Lymphoid$Post <- readRDS(file.path(genen_dir, "Innate_Lymphoid_Post_Filt_TOM.rds"))
TOM_sc$B_cell$Pre <- readRDS(file.path(genen_dir, "B_cell_Pre_Filt_TOM.rds"))
TOM_sc$B_cell$Post <- readRDS(file.path(genen_dir, "B_cell_Post_Filt_TOM.rds"))
TOM_sc$Myeloid_APC$Pre <- readRDS(file.path(genen_dir, "Myeloid_APC_Pre_Filt_TOM.rds"))
TOM_sc$Myeloid_APC$Post <- readRDS(file.path(genen_dir, "Myeloid_APC_Post_Filt_TOM.rds"))

length(intersect(colnames(TOM_sc$Myeloid_APC$Post), colnames(TOM_sc$Myeloid_APC$Pre)))

##########################################
### QUANTIFY REWIRING from pre to post ###
##########################################

analyze_rewiring <- function(TOM_pre_blocks, TOM_post_blocks, genes, top_k = 50) {
  res = list()
  for (gene in genes) {
    # figure out which gene is found in which block
    if (length(TOM_pre_blocks) == 2) {
      TOM_pre = find_gene_block(TOM_pre_blocks, gene)
      TOM_post = find_gene_block(TOM_post_blocks, gene)
    } else {
      TOM_pre = TOM_pre_blocks
      TOM_post = TOM_post_blocks
      #message(nrow(TOM_pre), " ", nrow(TOM_post))
    }
    
    # enforce shared indices
    tmp <- restrict_to_shared_genes(TOM_pre, TOM_post, gene)
    TOM_pre <- tmp$tom1
    TOM_post <- tmp$tom2
    
    # Get the top neighbors for each
    if (gene %in% intersect(rownames(TOM_pre), rownames(TOM_post))) {
      n_pre  <- get_top_neighbors(TOM_pre, gene, top_k = top_k)
      n_post <- get_top_neighbors(TOM_post, gene, top_k = top_k)
      
      # Jaccard overlap of neighbors
      jac <- jaccard(n_pre, n_post)
      
      # correlation of full row vectors (global similarity)
      row_cor <- safe_cor(TOM_pre[gene, ], TOM_post[gene, ], method = "spearman")
      
      # connectivity shift
      dk <- row_connectivity(TOM_post, gene) - row_connectivity(TOM_pre, gene)
    } else {
      jac = NA; row_cor = NA; dk = NA
    }
    
    
    res[[gene]] = data.frame(
      Gene = gene,
      jaccard_topK = jac,
      row_cor = row_cor,
      delta_connectivity = dk,
      stringsAsFactors = FALSE
    )
    
  }
  do.call(rbind, res)
}

bulk_rewiring <- analyze_rewiring(TOM_bulk$Pre, TOM_bulk$Post, genes_of_interest, top_k = 50)
bulk_rewiring$dataset <- "Bulk"
rownames(TOM_sc[["EOC"]]$Pre)
rownames(TOM_sc[["EOC"]]$Post)

sc_rewiring <- lapply(names(TOM_sc), function(ct) {
  df <- analyze_rewiring(TOM_sc[[ct]]$Pre, TOM_sc[[ct]]$Post, genes_of_interest, top_k = 50)
  df$dataset <- ct
  df
})
sc_rewiring <- do.call(rbind, sc_rewiring)

rewiring_all <- rbind(bulk_rewiring, sc_rewiring)

##################################################
### WHICH CELLTYPES EXPLAIN BULK REWIRING BEST ###
##################################################
compare_bulk_to_celltypes <- function(TOM_bulk_pre, TOM_bulk_post, TOM_sc_list, genes, sim_use="Corr") {
  out <- list(); out_pre <- list(); out_post <- list()
  for (gene in genes) {
    # figure out which gene is found in which block
    if (length(TOM_bulk_pre) == 2) {
      TOM_bulk_pre = find_gene_block(TOM_bulk_pre, gene)
      TOM_bulk_post = find_gene_block(TOM_bulk_post, gene)
    } else {
      TOM_bulk_pre = TOM_bulk_pre
      TOM_bulk_post = TOM_bulk_post
    }
   
    
    # for each cell type
    for (ct in names(TOM_sc_list)) {
      TOM_ct_pre <- TOM_sc_list[[ct]]$Pre
      TOM_ct_post <- TOM_sc_list[[ct]]$Post
      ###### DELTA ######
      g_shared_test <- Reduce(intersect, list(
        rownames(TOM_bulk_pre),
        rownames(TOM_bulk_post),
        rownames(TOM_ct_pre),
        rownames(TOM_ct_post)
      ))
      
      if (gene %in% g_shared_test) {
        # get teh shared gene universe between bulk and this cell type
        g_shared <- Reduce(intersect, list(
          colnames(TOM_bulk_pre),
          colnames(TOM_bulk_post),
          colnames(TOM_ct_pre),
          colnames(TOM_ct_post)
        ))
        cat("\nNumber DELTA shared genes for gene ", g, "\n in cell type ", ct, ": ", length(g_shared))
        # Get the change in TOM values for each gene (for universe)
        delta_bulk <- TOM_bulk_post[gene, g_shared] - TOM_bulk_pre[gene, g_shared]
        delta_ct   <- TOM_ct_post[gene, g_shared]   - TOM_ct_pre[gene, g_shared]
        # Get the correlation of the change
        if (sim_use == "jac") {
          sim <- NA
        } else {
          sim <- safe_cor(delta_bulk, delta_ct, method = "spearman")
        }
        
        out[[length(out) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = sim, stringsAsFactors = FALSE
        )
      } else {
        out[[length(out) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = NA, stringsAsFactors = FALSE
        )
      }
      ###### PRE ######
      if (gene %in% intersect(rownames(TOM_bulk_pre), rownames(TOM_ct_pre))) {
        # restrict to shared genes space
        tmp <- restrict_to_shared_genes(TOM_bulk_pre, TOM_ct_pre, gene)
        TOM_bulk_pre_tmp <- tmp$tom1
        TOM_ct_pre_tmp <- tmp$tom2
        if (sim_use == "jac") {
          n_bulk  <- get_top_neighbors(TOM_bulk_pre_tmp, gene, top_k = 100)
          n_ct <- get_top_neighbors(TOM_ct_pre_tmp, gene, top_k = 100)
          # Jaccard overlap of neighbors
          sim <- jaccard(n_bulk, n_ct)
        } else {
          # Get the correlation of the change
          sim <- safe_cor(TOM_bulk_pre_tmp[gene,], TOM_ct_pre_tmp[gene,], method = "spearman")
        }
        out_pre[[length(out_pre) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = sim, stringsAsFactors = FALSE
        )
      } else {
        out_pre[[length(out_pre) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = NA, stringsAsFactors = FALSE
        )
      }
      ###### POST ######
      if (gene %in% intersect(rownames(TOM_bulk_post), rownames(TOM_ct_post))) {
        # restrict to shared genes space
        tmp <- restrict_to_shared_genes(TOM_bulk_post, TOM_ct_post, gene)
        TOM_bulk_post_tmp <- tmp$tom1
        TOM_ct_post_tmp <- tmp$tom2
        if (sim_use == "jac") {
          n_bulk  <- get_top_neighbors(TOM_bulk_post_tmp, gene, top_k = 100)
          n_ct <- get_top_neighbors(TOM_ct_post_tmp, gene, top_k = 100)
          # Jaccard overlap of neighbors
          sim <- jaccard(n_bulk, n_ct)
        } else {
          # Get the correlation of the change
          sim <- safe_cor(TOM_bulk_post_tmp[gene,], TOM_ct_post_tmp[gene,], method = "spearman")
        }
        out_post[[length(out_post) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = sim, stringsAsFactors = FALSE
        )
      } else {
        out_post[[length(out_post) + 1]] <- data.frame(
          Gene = gene, CellType = ct,
          TOM_similarity = NA, stringsAsFactors = FALSE
        )
      }
    }
  }
  out = do.call(rbind, out)
  out_pre = do.call(rbind, out_pre)
  out_post = do.call(rbind, out_post)
  return(list("Delta"=out, "Pre"=out_pre, "Post"=out_post))
  
}

bulk_celltype_similarity <- compare_bulk_to_celltypes(
  TOM_bulk$Pre, TOM_bulk$Post,
  TOM_sc,
  genes_of_interest
)

bulk_celltype_jaccard <- compare_bulk_to_celltypes(
  TOM_bulk$Pre, TOM_bulk$Post,
  TOM_sc,
  genes_of_interest, sim_use="jac"
)

get_celltype_summary <- function(data_use) {
  celltype_summary <- aggregate(
    TOM_similarity ~ CellType,
    data = data_use,
    FUN = function(x) c(mean = mean(x, na.rm=TRUE), median = median(x, na.rm=TRUE), n = sum(!is.na(x)))
  )
  
  # unpack matrix columns
  celltype_summary <- data.frame(
    CellType = celltype_summary$CellType,
    mean_similarity = celltype_summary$TOM_similarity[, "mean"],
    median_similarity = celltype_summary$TOM_similarity[, "median"],
    n_genes = celltype_summary$TOM_similarity[, "n"]
  )
  
  celltype_summary <- celltype_summary[order(celltype_summary$mean_similarity, decreasing = TRUE), ]
  return(celltype_summary)
  
}

CT_delta_summary = get_celltype_summary(bulk_celltype_similarity$Delta)
CT_pre_summary = get_celltype_summary(bulk_celltype_similarity$Pre)
CT_post_summary = get_celltype_summary(bulk_celltype_similarity$Post)

################
### PLOTTING ###
################

# Heatmap for CT vs Bulk Similarity
library(ComplexHeatmap)
library(circlize)

plot_Heatmap <- function(B_CT_sim, title_use) {
  
  mat <- reshape(
    B_CT_sim,
    idvar = "Gene",
    timevar = "CellType",
    direction = "wide"
  )
  
  rownames(mat) <- mat$Gene
  mat$Gene <- NULL
  mat <- as.matrix(mat)
  
  # clean colnames
  colnames(mat) <- gsub("^TOM_similarity\\.", "", colnames(mat))
  
  # -------------------------
  # 1) REMOVE ROWS WITH ALL NA
  # -------------------------
  keep_rows <- rowSums(!is.na(mat)) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  
  if (nrow(mat) == 0) {
    stop("After removing all-NA rows, no genes remain to plot.")
  }
  
  # -------------------------
  # 2) LABEL ONLY THE MAX CELL PER ROW
  # -------------------------
  # create a label matrix (same dim as mat)
  label_mat <- matrix(
    "",
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  # for each row, find the max (ignoring NA), label only that cell
  for (i in seq_len(nrow(mat))) {
    row_vals <- mat[i, ]
    if (all(is.na(row_vals))) next
                                  
    j_max <- which.max(row_vals)  # which.max ignores NA by treating as -Inf effectively
    label_mat[i, j_max] <- sprintf("%.2f", row_vals[j_max])
    # get the values above 0.1 too
    j_abov1 <- which(row_vals > 0.1)
    label_mat[i, j_abov1] <- sprintf("%.2f", row_vals[j_abov1])
  }
  
  if (min(mat, na.rm=TRUE) < 0) {
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  } else {
    col_fun <- colorRamp2(c(0, 0.4, 0.8), c("blue", "white", "red"))
  }
  
  
  Heatmap(
    mat,
    name = title_use,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "grey85"),
    column_names_rot = 90,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- label_mat[i, j]
      if (lab != "") {
        grid.text(
          lab,
          x, y,
          gp = gpar(fontsize = 10, col = "black", fontface = "bold")
        )
      }
    }
  )
}


plot_Heatmap(bulk_celltype_similarity$Delta, "∆TOM Correlation\n(bulk vs celltype)")
plot_Heatmap(bulk_celltype_similarity$Pre, "Pre TOM Correlation\n(bulk vs celltype)")
plot_Heatmap(bulk_celltype_similarity$Post, "Post TOM Correlation\n(bulk vs celltype)")

plot_Heatmap(bulk_celltype_jaccard$Pre, "Pre TOM\nTop 100 Jaccard\n(bulk vs celltype)")
plot_Heatmap(bulk_celltype_jaccard$Post, "Post TOM\nTop 50 Jaccard\n(bulk vs celltype)")



