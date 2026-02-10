## Considering Longitudinal/Non-Longitudinal Microarray & scRNA-seq data for HGSOC

### Prepping data for analysis (SC_Micro)
1. Prep the scRNA-seq (and patient-matched Bulk RNA-seq) for downstream comparison in `00_Prep_scRNAseq.ipynb`

### Analysis
#### scRNA-seq (SC_Micro)
1. Compare patient-matched scRNA-seq and RNA-seq in both expression and log2FC in `01_Compare_RNA_scRNAseq.ipynb`
2. Assess cell-type specific expression and association, and if cell type proportions are associated with PFI in `02_RNA_Celltype.ipynb`

#### Non-longitudinal Microarray (SC_Micro)
1. Get non-longitudinal microarray data and perform univariate Cox regression analyses & Kaplan-Meier curves in `03_Format_Microarray_KM_curves.ipynb` :: NOTE: need to clean up since LOTS of plots

#### Network analysis (Network)
1. `04_hdWGCNA_Networks.ipynb`: Prepare the Bulk RNA-seq for WGCNA, run hdWGCNA on single-cell data and save TOMs and annotate modules.
2. `05_WGCNA_Networks.R`: Run WGCNA on Bulk RNA-seq and save TOMs and annotate modules.
3. `06_Network_Compare.R`: Compare the networks of the model genes in bulk to those in scRNA-seq (at cell type level)
4. `07_GBP4_Networks.ipynb`: Identifying and plotting GBP4 networks based on fucntional annotation and timing-specificity.