# Make sure it is considering the correct libraries)
#.libPaths(c("/Users/hoto7260/R/GLMNET", "/Users/hoto7260/R/x86_64-pc-linux-gnu-library/4.4", "/opt/R/4.4.0/lib/R/library/"))
.libPaths(c("/Users/hoto7260/R/GLMNET",  "/opt/R/4.4.0/lib/R/library/"))


library(WGCNA)



#library(data.table)
#library(ggplot2)
#library(grDevices)
#library(DESeq2)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
#install.packages("igraph")
library(igraph)
#install.packages("BiocManager") # vs 3.20
# BiocManager::install("fgsea") <- had to allow it to be installed in my 4/x86_64-pc one despite it using the same GLMNET
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#library(clusterProfiler)
library(WGCNA)
library(fgsea, lib.loc=c("/Users/hoto7260/R/x86_64-pc-linux-gnu-library/4.4"))




allowWGCNAThreads(24)  # Enable multi-threading

# Set working directory or output path as needed
plot_dir <- "/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/plots"
dir.create(plot_dir, showWarnings = FALSE)

gene_list = c("IDO1", "CD40", "JAK2", "CTSS", "ALDOC", "TMUB2", "ID4", "ENO1", "EGR1", "OAS3", 
              "RRM2", "GBP4", "AMOTL2", "ESYT3", "HCN3", "DDX11", "GPR173", "CLK2", "CDKL2", "SHROOM1", "CD274")
length(gene_list)


######################################################################
#### 1. Load and merge RNA-seq datasets ####
######################################################################
timing = "POST"
datExpr_bulk <- readRDS(paste0("/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/counts/WGCNA_", timing, "_Bulk_VST_data.rds"))
dim(datExpr_bulk)


spt <- pickSoftThreshold(datExpr_bulk) 

bulk_net <- blockwiseModules(
  datExpr_bulk,
  power = spt$powerEstimate,
  loadTOM=TRUE,
  maxBlockSize = 25000,
  TOMType = "signed",
  minModuleSize = 30,
  randomSeed = 20,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE, # saves TOM matrices to disk to avoid recomputing
  saveTOMFileBase = paste0("TOM_Full_", timing, "_block"), # prefix for saved files
  verbose = 3
)

bulk_colors <- labels2colors(bulk_net$colors)
bulk_MEs <- bulk_net$MEs

######################################################################
#### 3. Evaluate Module members ####
######################################################################

bulk_MM <- cor(datExpr_bulk, bulk_MEs, use = "p")
bulk_MM_p <- corPvalueStudent(bulk_MM, nSamples = nrow(datExpr_bulk))


# Subset and coerce to data.frame
bulk_mm_sub <- as.data.frame(bulk_MM[gene_list, , drop = FALSE])
# Add Gene column
bulk_mm_sub$Gene <- rownames(bulk_mm_sub)
# Convert to long format manually (tidyverse isn't installing)
bulk_mm_long <- data.frame(
  Gene   = rep(bulk_mm_sub$Gene, times = ncol(bulk_mm_sub) - 1),
  Module = rep(colnames(bulk_mm_sub)[colnames(bulk_mm_sub) != "Gene"],
               each = nrow(bulk_mm_sub)),
  MM     = as.vector(as.matrix(bulk_mm_sub[, colnames(bulk_mm_sub) != "Gene"]))
)

# Add p-values (same order as MM vector)
bulk_mm_long$pval <- as.vector(bulk_MM_p[gene_list, , drop = FALSE])
# Add source column
bulk_mm_long$source <- "Bulk"

mm_mat <- with(
  bulk_mm_long,
  tapply(MM, list(Gene, Module), identity)
)

mm_mat <- as.matrix(mm_mat)


Heatmap(
  mm_mat,
  col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
  name = "Correlation to Modules",
  row_title = "Model Genes",
  column_title = "WGCNA networks",
  cluster_rows = TRUE,
  cluster_columns = TRUE
)


# graph again but only considering modules with the gene list

# Extract module assignments
module_colors <- bulk_net$colors
names(module_colors) <- colnames(datExpr_bulk)
module_colors[gene_list]
gene_modules <- unique(module_colors[gene_list])
length(gene_list)
length(gene_modules)
length(unique(module_colors))

mm_mat <- with(
  bulk_mm_long[bulk_mm_long$Module %in% paste0("ME", gene_modules),],
  tapply(MM, list(Gene, Module), identity)
)


mm_mat <- as.matrix(mm_mat)
pdf(paste0("Cell_Module_", timing, "_Heatmap.pdf"), width = 7, height = 7)

ht = Heatmap(
  mm_mat,
  col = colorRamp2(c(0, 0.5, 1), c("#417ca8ff", "white", "#da3a46ff")),
  name = "Module\ncorrelation",
  row_title = "Model Genes",
  column_title = "WGCNA networks",
  cluster_rows = TRUE,
  cluster_columns = TRUE, 
  cell_fun = function(j, i, x, y, width, height, fill) {
    gene <- rownames(mm_mat)[i]
    module <- colnames(mm_mat)[j]
    gene_module = paste0("ME",module_colors[gene])
    if (!is.na(gene_module) &&
        gene_module == module) {
      
      grid.text(
        "X",
        x = x,
        y = y,
        gp = gpar(col = "black", fontsize = 14, fontface = "bold")
      )
    }
  }
)
print(ht)
dev.off()
######################################################################
#### 4. Biologically annotate each module with a core gene in it ####
######################################################################
# Get the gene lists per module
module_genes <- lapply(gene_modules, function(m) {
  names(module_colors)[module_colors == m]
})
names(module_genes) <- gene_modules

### Perform GO enrichment
#BiocManager::install("fgsea", version="3.20")
annotate_module_GO <- function(genes, universe) {
  
  message(length(genes))
  
  eg <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  eg <- na.omit(eg)
  
  universe_eg <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = universe,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  universe_eg <- na.omit(universe_eg)
  
  clusterProfiler::enrichGO(
    gene = eg,
    universe = universe_eg,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
}

# Annotate all modules
module_GO <- lapply(
  module_genes, annotate_module_GO, universe = colnames(datExpr_bulk)
)

# Print out the top processes per module
module_summary <- lapply(module_GO, function(x) {
  if (is.null(x) || nrow(x@result) == 0) return(NA)
  head(x@result[, c("Description", "p.adjust", "geneID")], 10)
})
module_summary

full_module_go_names <- list("30"="lipid|Carb_transport", "13"="cell_fate", 
                        "2"="expr_regulation", "15"="Unclear", "3"="Macroautophagy", 
                        "7"="Mitosis", "4"="neuron|muscle|GProtein", "1"="Mitochon", 
                        "5"="RNA_splicing", "0"="chemic_stimulus|bact_response", 
                        "9"="lymph/leuk_activation", "8"="leuk_med_killing", "45"="viral_response")

full_module_go_names <- list("30"="TollLR_Leuk", "13"="sex_diff", 
                        "2"="expr_regulation", "15"="ER_transport", "3"="Prot_local", 
                        "7"="Mitosis", "4"="Muscle", "1"="RNAsplic|Mitochon", 
                        "5"="ER_stress|ERAD", "0"="chemic_stimulus|Metab", 
                        "9"="lymph/leuk_activation", "8"="leuk_med_killing", "45"="viral_response")

pre_module_go_names <- list("20"="neuron_diff", "3"="Respiration", "5"="Nucleosome_DNA", 
                            "2"="macroautophagy_cellpolarity", "15"="Mitosis", "1"="Splicing_Degrad", 
                            "13"="apoptosis_signal", "0"="chemstimulus_humoralimmune", 
                            "9"="ImmuneReceptor_LeukAdhesion", "10"="LeukMig_Activation", 
                            "45"="viral_response")

post_module_go_names <- list("23"="plasma_lipoprotein", "10"="Mix_protein_ER", "2"="macroautophagy_Transcrip", 
                             "5"="MitRespir_DNArepair", "12"="CellCycle", "1"="ChemStim_mRNAdegrad",
                             "3"="MHC_ER", "13"="EpithelCytoskeletondvlpt", "16"="ImmuneRecpt_TNF", 
                             "4"="Cilia_Muscle", "9"="AdaptiveImmune")

module_colors[gene_list]

module_summary <- lapply(module_GO, function(x) {
  if (is.null(x) || nrow(x@result) == 0) return(NA)
  x_filt <- x@result[grepl(paste(gene_list, collapse="|"), x@result$geneID),]
  head(x_filt[, c("Description", "p.adjust", "geneID")], 10)
})
module_summary


######################################################################
#### 5. Explore gene networks of model genes for cytoscape
###################################################################
# create output directory
genen_dir <- "/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/gene_networks"
dir.create(genen_dir, showWarnings = FALSE)

build_gene_network_from_TOM <- function(
    gene,
    TOM_sim,
    top_n = 30,
    tom_threshold = 0.8
) {
  if (!gene %in% colnames(TOM_sim)) return(NULL)
  
  tom_vec <- TOM_sim[gene, ]
  tom_vec <- tom_vec[names(tom_vec) != gene]
  # only keep the non pseudogenes, ribosomal, mitochondrial, lncRNA
  tom_vec <- tom_vec[!grepl("^(LIN|LOC|SNOR|MT-|RPL|RPS|RP)", names(tom_vec))]
  max_sim = max(tom_vec)
  message(paste0("Max similarity:", max_sim))
  
  keep <- which(tom_vec >= tom_threshold)
  if (length(keep) == 0) return(NULL)
  
  tom_vec <- sort(tom_vec[keep], decreasing = TRUE)
  
  #message(paste(gene, tom_vec["GBP5"], tom_vec["GBP4"], tom_vec["GBP1"]))
  
  message(paste0("\nNumber of genes with similarity >", tom_threshold, ": ", length(keep)))
  if (length(tom_vec) > top_n) {
    tom_vec <- tom_vec[seq_len(top_n)]
  }
  
  edge_df <- data.frame(
    from = gene,
    to = names(tom_vec),
    weight = as.numeric(tom_vec),
    stringsAsFactors = FALSE
  )
  edge_df$Type <- ">0.05"
  if (max_sim > 0.1) {
    edge_df[edge_df$weight > 0.1,]$Type <- ">0.1"
  }
  if (max_sim > 0.2) {
    edge_df[edge_df$weight > 0.2,]$Type <- ">0.2"
  }
  if (max_sim > 0.25) {
    edge_df[edge_df$weight > 0.25,]$Type <- ">0.25"
  }
  
  
  node_df <- unique(
    rbind(
      data.frame(name = gene, type = "ModelGene"),
      data.frame(name = edge_df$to, type = "Neighbor")
    )
  )
  
  list(edges = edge_df, nodes = node_df)
}

gene_net_go <- list()
full_gene_list <- c()
TOM_keep_list <- list()

for (b in c("1", "2")) {
  message("Loading block ", b)
  load(paste0("TOM_Full_", timing, "_block-block.", b, ".RData"))  # loads TOM
  
  # convert TOM to matrix
  # Map indices → gene names
  block_idx   <- bulk_net$blockGenes[[as.integer(b)]]
  block_genes <- colnames(datExpr_bulk)[block_idx]
  
  # Convert TOM from dist → full matrix
  TOM_dist <- as.dist(TOM)
  TOM_mat <- as.matrix(TOM)
  colnames(TOM_mat) <- block_genes
  rownames(TOM_mat) <- block_genes
  message(block_genes[1:4])
  message(colnames(TOM)[1:4])
  
  
  
  genes_in_block <- intersect(
    gene_list,
    colnames(TOM_mat)
  )
  
  for (g in genes_in_block) {
    
    message("Building network for ", g)
    net <- build_gene_network_from_TOM(
      gene = g,
      TOM = TOM_mat,
      top_n = 30,
      tom_threshold = 0.05
    )
    
    if (is.null(net)) next
    
    
    # ---- Save Cytoscape files ----
    write.table(
      net$edges,
      file = file.path(genen_dir, paste0(g, "_", timing, "_edges.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    # print out the summary GO
    #gene_net_go[g] = annotate_module_GO(c(g, net$edges$to), block_genes)
    # add to gene list 
    full_gene_list = c(full_gene_list, g, net$edges$to)
    # remove any difficult to characterize genes
    full_gene_list <- full_gene_list[!grepl("^(LIN|LOC|SNOR|MT-|RPL|RPS|RP)", full_gene_list)]
    
  }
  # Get the TOM parts for the genes of interest
  genes_consider = intersect(colnames(TOM_mat), full_gene_list)
  TOM_keep_list[[b]] <- TOM_mat[genes_consider,genes_consider]
}

# get only the genes considered in more than one network
dup_full_gene_list <- union(full_gene_list[duplicated(full_gene_list)], gene_list)
length(dup_full_gene_list)
# get all genes
full_gene_list <- unique(full_gene_list)
length(full_gene_list)

gene_functional_annotation <- list("T cell activation"=c("CD274", "SLAMF6", "CCR2", "THEMIS", "TNFSF8","CD3G","TIGIT","JAML","KLRK1","KLRC4-KLRK1","SASH3"),
                                   "Type II IFN"=c("CD274","SLAMF6","CCR2","KLRK1","KLRC4-KLRK1","SASH3","GBP4","CCL5"), 
                                   "Response to nitrogen or hormone"=c("GPR173", "SHOC2","OSBPL8","RAP1A","AGRN"), 
                                   "Cell cycle"=c("RRM2","NUF2","KIF2C","KIF14","BUB1","KIF11","MELK","NEK2","RACGAP1","DLGAP5","CCNB2","CEP55","CCNB1",
                                                  "KNL1","KIF4A","BUB1B","KIF20A","NCAPG","KIF15"), 
                                   "Immune Response"=c("GBP4","IL10RA","MEPEG1","CD96","GZMA","OAS3","SASH3","LILRB1","CD4","LILRB4","NCKAP1L","HAVCR2","LAPTM5",
                                                       "CTSS","CCL5","CD86","ITGB2","SLAMF6"), 
                                   "TNF cytokine prod"=c("OAS3","SASH3","LILRB1","LILRB4","HAVCR2","CD86",'OAS2'), 
                                   "muscle structure"=c("EGR1","SYNE1","LAMA2","MYH11","LMOD1","DMD"), 
                                   "muscle contraction"=c("ENO1","MYH11","LMOD1","DMD"), 
                                   "proteasome"=c("TMUB2","TMEM259","SIRT6","FBXL15"), 
                                   "Lymph/Leuk activation"=c("CD40","ITGAL","SLAMF6","SASH3","DOCK2","CD96","LY9","PIK3CG","SLAMF1","CD2","TRAT1","CD3D","ITK","CCL5"), 
                                   "Leuk cell-cell adhesion"=c("IDO1","CCR2","TIGIT","LAX1","SASH3","CCL5","ICOS"))
# Pre
gene_functional_annotation <- list(
  "Immune Response"=c("CTSS", "CD300A", "HAVCR2", "TLR2", "TLR4", "PTPRC", "AIM2", "SLAMF6", 
                      "TNFSF13B", "CLEC7A", "TLR1", "LY96", "CD86", "C3AR1", "CCR2", "CASP1", "PIK3CG", 
                      "OAS3", "KLHL6", "LILRB1", "SLAMF6", "CCR2", "SASH3", "CCL5", "CD38", "LILRB4", 
                      "LAPTM5", "CD28", "FPR1", "TLR2", "FCER1G", "CD33", "CLEC7A", "FCGR3A", 
                      "LAIR1", "VSIG4","CRTAM", "IDO1", "AIM2", "TNFSF13B", "CASP1"), 
  "Lymph/Leuk Activation"=c("CD274", "TNFSF13B", "PDCD1LG2", "SLAMF6", "CD86", "TLR4", "CD84", 
                            "IL2RA", "TNFSF8", "IL7R", "TIGIT", "CLEC7A", "FGL2", "JAK2", "PIK3CG", 
                            "PTPN22", "PTPRC", "TLR4", "FGL2", "DOCK2", "BTK", "SAMSN1", "NCKAP1L", 
                            "CD40", "DOCK2", "PIK3CG", "NCKAP1L", "ITK", "SLAMF6", 
                            "ITGAL", "LILRB1", "JAML", "CD2", "PRKCB", "IKZF1"), 
  "Lymph/Leuk/Mono/Mye proliferation"=c("JAK2", "PIK3CG", "PTPN22", "PTPRC", "TLR4", "DOCK2", "BTK", "NCKAP1L", 
                                        "LCP2", "CD40", "LILRB1", "CSF2RB"), 
  "DNA damage"=c("DDX11", "POLD1", "BRAT1", "KMT5C", "TELO2", "UPF1", "PIDD1", "ATAD3A", "XRCC3"), 
  "Mitosis"=c("RRM2", "NUF2", "KIF11", "DBF4", "CCNA2", "RACGAP1", "CDK1", "CDKN3", "TTK", "CCNB1", 
              "CKAP2", "NCAPG", "NEK2", "NUSAP1", "CEP55", "KIF15", "DLGAP5", "VRK1"), 
  "TNF cytokine"=c("OAS3", "IL10RA", "LILRB1", "CCR2", "CSF2RB", "CCL5", "CCR5", "IL21R", 
                   "LILRB4", "SASH3", "GPR174", "TIGIT", "SLAMF6", "JAK2", "PTPN22", 
                   "PTPRC", "TLR4", "LY96", "BTK"), 
  "nucleobase-catabolism"=c("ENO1", "ALDOC", "NORAD", "OIP5-AS1")
)
# Post
gene_functional_annotation <- list(
  "Immune Process"=c("CD274", "CD226", "IRAK3", "CD84", "KLRD1", "SH2D1A", "TLR4", "SLAMF8", "FCGR2B", 
                     "GBP4", "IL10RA", "KLRC4-KLRK1", "CCR5", "KLRK1", "SPN", "MPEG1", "GZMA", 
                     "OAS3", "CD4", "SPN", "LILRB1", "CSF2RB", "KLHL6", "CTSS", "LAIR1", "NCKAP1L", 
                     "HAVCR2", "CD300A", "LILRB4", "CLEC7A", "CD300LF", "CD40", "SLAMF6", 
                     "SH2D1A", "SASH3", "SLAMF1", "PRF1"), 
  "Lymph/Leuk Activation"=c("CD274", "TNFSF13B", "TNFSF8", "CD226", "CD84", "KLRD1", "TLR4", "SLAMF8", "FCGR2B", 
                            "CD40", "CCL5", "SASH3", "SLAMF1", "SPN", "IDO1", "SLAMF6", "IKZF3", "CCR2", "NCKAP1L", 
                            "CD3G", "CD86", "SLAMF1", "CD3D", "TIGIT", "JAML", "CCL5"), 
  "Lymph/Leuk/Mono/Mye proliferation"=c("IDO1", "SASH3", "IKZF3", "CCR2", "SPN", "NCKAP1L", "CD86", "SLAMF1", "CCL5", 
                                        "CD40", "CCL5", "SASH3", "SLAMF1", "SPN"), 
  "DNA damage"=c("DDX11", "POLD1", "TELO2", "SIRT6", "AP5Z1"), 
  "Mitosis"=c("RRM2", "ECT2", "NUF2", "BUB1", "KIF11", "NEK2", "KIF14", "KIF2C", "MELK", "TTK", "CEP55", "CLSPN", 
              "KIF4A", "CENPF", "BRCA1", "BUB1B", "DLGAP5", "CCNB2", "NCAPH"), 
  "TNF cytokine"=c("OAS3", "IL10RA", "CD4", "LILRB1", "IL12RB1", "CSF2RB", "CCL5", "SASH3", "SPN"), 
  "Epithelial differentiation"=c("ALDOC", "OVOL2", "ESRP1", "EVPL", "CLDN3", "GRHL2", "KDF1", "PAX8", "TJP3")
  
)

gene_functional_annotation <- list(
  "Immune Process"=c("CD274", "CD226", "IRAK3", "CD84", "KLRD1", "SH2D1A", "TLR4", "SLAMF8", "FCGR2B", 
                     "GBP4", "IL10RA", "KLRC4-KLRK1", "CCR5", "KLRK1", "SPN", "MPEG1", "GZMA", 
                     "OAS3", "CD4", "SPN", "LILRB1", "CSF2RB", "KLHL6", "CTSS", "LAIR1", "NCKAP1L", 
                     "HAVCR2", "CD300A", "LILRB4", "CLEC7A", "CD300LF", "CD40", "SLAMF6", 
                     "SH2D1A", "SASH3", "SLAMF1", "PRF1", 
                     "CTSS", "CD300A", "HAVCR2", "TLR2", "TLR4", "PTPRC", "AIM2", "SLAMF6", 
                     "TNFSF13B", "CLEC7A", "TLR1", "LY96", "CD86", "C3AR1", "CCR2", "CASP1", "PIK3CG", 
                     "OAS3", "KLHL6", "LILRB1", "SLAMF6", "CCR2", "SASH3", "CCL5", "CD38", "LILRB4", 
                     "LAPTM5", "CD28", "FPR1", "TLR2", "FCER1G", "CD33", "CLEC7A", "FCGR3A", 
                     "LAIR1", "VSIG4","CRTAM", "IDO1", "AIM2", "TNFSF13B", "CASP1"), 
  "Lymph/Leuk Activation"=c("CD274", "TNFSF13B", "TNFSF8", "CD226", "CD84", "KLRD1", "TLR4", "SLAMF8", "FCGR2B", 
                            "CD40", "CCL5", "SASH3", "SLAMF1", "SPN", "IDO1", "SLAMF6", "IKZF3", "CCR2", "NCKAP1L", 
                            "CD3G", "CD86", "SLAMF1", "CD3D", "TIGIT", "JAML", "CCL5", "CD274", "TNFSF13B", "PDCD1LG2", "SLAMF6", "CD86", "TLR4", "CD84", 
                            "IL2RA", "TNFSF8", "IL7R", "TIGIT", "CLEC7A", "FGL2", "JAK2", "PIK3CG", 
                            "PTPN22", "PTPRC", "TLR4", "FGL2", "DOCK2", "BTK", "SAMSN1", "NCKAP1L", 
                            "CD40", "DOCK2", "PIK3CG", "NCKAP1L", "ITK", "SLAMF6", 
                            "ITGAL", "LILRB1", "JAML", "CD2", "PRKCB", "IKZF1"), 
  "Lymph/Leuk/Mono/Mye proliferation"=c("IDO1", "SASH3", "IKZF3", "CCR2", "SPN", "NCKAP1L", "CD86", "SLAMF1", "CCL5", 
                                        "CD40", "CCL5", "SASH3", "SLAMF1", "SPN", "JAK2", "PIK3CG", "PTPN22", "PTPRC", "TLR4", "DOCK2", "BTK", "NCKAP1L", 
                                        "LCP2", "CD40", "LILRB1", "CSF2RB"), 
  "DNA damage"=c("DDX11", "POLD1", "TELO2", "SIRT6", "AP5Z1", "DDX11", "POLD1", "BRAT1", "KMT5C", "TELO2", "UPF1", "PIDD1", "ATAD3A", "XRCC3"), 
  "Mitosis"=c("RRM2", "ECT2", "NUF2", "BUB1", "KIF11", "NEK2", "KIF14", "KIF2C", "MELK", "TTK", "CEP55", "CLSPN", 
              "KIF4A", "CENPF", "BRCA1", "BUB1B", "DLGAP5", "CCNB2", "NCAPH", "RRM2", "NUF2", "KIF11", "DBF4", "CCNA2", "RACGAP1", "CDK1", "CDKN3", "TTK", "CCNB1", 
              "CKAP2", "NCAPG", "NEK2", "NUSAP1", "CEP55", "KIF15", "DLGAP5", "VRK1"), 
  "TNF cytokine"=c("OAS3", "IL10RA", "CD4", "LILRB1", "IL12RB1", "CSF2RB", "CCL5", "SASH3", "SPN", 
                   "OAS3", "IL10RA", "LILRB1", "CCR2", "CSF2RB", "CCL5", "CCR5", "IL21R", 
                   "LILRB4", "SASH3", "GPR174", "TIGIT", "SLAMF6", "JAK2", "PTPN22", 
                   "PTPRC", "TLR4", "LY96", "BTK"), 
  "Epithelial differentiation"=c("ALDOC", "OVOL2", "ESRP1", "EVPL", "CLDN3", "GRHL2", "KDF1", "PAX8", "TJP3"), 
  "nucleobase-catabolism"=c("ENO1", "ALDOC", "NORAD", "OIP5-AS1")
  
)

edge_df <- data.table::fread("/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/gene_networks/POST_FULL_edges.tsv")
edge_df[1:2,]
dim(edge_df)
edge_df <- edge_df[(edge_df$From %in% gene_list) | (edge_df$To %in% gene_list),]
dim(edge_df)
edge_df[1:2,]
write.table(edge_df, "/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/gene_networks/POST_FULL_FILT_edges.tsv", row.names=FALSE, quote=FALSE, sep="\t")
node_df <- data.table::fread("~/POST_FULL_nodes.tsv", fill=TRUE)
node_df[1:2,]
node_df$Model_Gene <- FALSE
node_df[node_df$Gene %in% gene_list,]$Model_Gene <- TRUE

cat_gene_list = c()
for (gene in node_df$Gene) {
  category = "None"
  # get the category
  for (annotation in names(gene_functional_annotation)) {
    #print(gene_functional_annotation[[annotation]])
    if (gene %in% gene_functional_annotation[[annotation]]) {
      #message(paste(gene, "in", annotation))
      if (category == "None") {
        category = annotation
      } else {
        category = paste0(category, "|", annotation)
      }
    }
  }
  cat_gene_list = c(cat_gene_list, category)
  
}
node_df$Category <- cat_gene_list
node_df[1:2,]
table(node_df$Category)
write.table(node_df, "/scratch/Shares/clauset/Clauset_ABNexus/WGCNA/gene_networks/POST_FULL_nodes.tsv", row.names=FALSE, quote=FALSE, sep="\t")









gene_module_summary <- lapply(gene_net_go, function(x) {
  if (is.null(x) || nrow(x@result) == 0) return(NA)
  x_filt <- x@result[grepl(paste(gene_list, collapse="|"), x@result$geneID),]
  head(x_filt[, c("Description", "p.adjust", "geneID")], 10)
})
gene_module_summary


intersect(rownames(TOM_keep_list[["1"]]), colnames(TOM_keep_list[["1"]]))
ncol(TOM_keep_list[["1"]])
nrow(TOM_keep_list[["1"]])

###################################
# Now create a full gene network (get all the edges)
TOM_use = TOM_keep_list[["1"]]
TOM_use[lower.tri(TOM_use)] <- 0 # Keep upper triangle only
diag(TOM_use) <- 0 # Remove self-loops
edgeList <- reshape2::melt(TOM_use)
colnames(edgeList) <- c("From", "To", "Weight")
edgeList <- edgeList[edgeList$Weight > 0.1, ] # Filter

TOM_use = TOM_keep_list[["2"]]
TOM_use[lower.tri(TOM_use)] <- 0 # Keep upper triangle only
diag(TOM_use) <- 0 # Remove self-loops
edgeList2 <- reshape2::melt(TOM_use)
colnames(edgeList2) <- c("From", "To", "Weight")
edgeList2 <- edgeList2[edgeList2$Weight > 0.1, ] # Filter

dim(edgeList)
dim(edgeList2)
edgeList <- rbind(edgeList, edgeList2)
dim(edgeList)
# remove any edges that aren't to/from a Model gene
edgeList <- edgeList[edgeList$From %in% gene_list | edgeList$To %in% gene_list,]
dim(edgeList)
write.table(
  edgeList,
  file = file.path(genen_dir, paste0(timing, "_FULL_edges.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



# add the categories of each gene -- if multiple then add
cat_gene_list = c()
for (gene in full_gene_list) {
  category = ""
  # get the category
  for (annotation in names(gene_functional_annotation)) {
    #print(gene_functional_annotation[[annotation]])
    if (gene %in% gene_functional_annotation[[annotation]]) {
      #message(paste(gene, "in", annotation))
      if (category == "") {
        category = annotation
      } else {
        category = paste0(category, "|", annotation)
      }
    }
  }
  cat_gene_list = c(cat_gene_list, category)
  
}

cat_gene_df = as.data.frame(data.table::data.table("Gene"=full_gene_list, "Category"=cat_gene_list))
cat_gene_df[cat_gene_df$Gene %in% gene_list,]

write.table(
  cat_gene_df,
  file = file.path(genen_dir, paste0(timing, "_FULL_nodes.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# only save the TOM for the gene list to correlate with the sc data
for (b in c("1", "2")) {
  message("Loading block ", b)
  load(paste0("TOM_Full_", timing, "_block-block.", b, ".RData"))  # loads TOM
  
  # convert TOM to matrix
  # Map indices → gene names
  block_idx   <- bulk_net$blockGenes[[as.integer(b)]]
  block_genes <- colnames(datExpr_bulk)[block_idx]
  
  # Convert TOM from dist → full matrix
  TOM_dist <- as.dist(TOM)
  TOM_mat <- as.matrix(TOM)
  colnames(TOM_mat) <- block_genes
  rownames(TOM_mat) <- block_genes
  message(block_genes[1:4])
  message(colnames(TOM)[1:4])
  
  
  genes_in_block <- intersect(
    gene_list,
    colnames(TOM_mat)
  )
  
  TOM_consider = TOM_mat[genes_in_block,]
  message(nrow(TOM_consider), " ", ncol(TOM_consider))
  saveRDS(TOM_consider, file.path(genen_dir, paste0("TOM_filt", b, "_", timing, ".tsv")))
  
}
TOM_keep_list[["1"]][intersect(rownames(TOM_keep_list[["1"]]), gene_list),]







