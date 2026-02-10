## Paired HGSOC Pre & Post NACT Sample Analysis

### Organizing the Metadata and NanoString Counts (Preprocessing)

Used to get NonEGA, and NanoString formatted Counts

1. Get all metadata into one file and format NanoString counts. Nanostring/analysis/`00_Get_metadata_and_format_counts.ipynb`
2. EGA data had several sources of metadata that needed to be linked together. This was done using `01_Harmonizing_Metadata.ipynb` --> should be EGA 
3. Combine all metadata with `02_Combine_EGA_Rest_Metadata.ipynb`

### Get RNA-seq Counts (Preprocessing)
Javella et al. only provided already normalized counts so we could not map and count for ourselves.
1. Get the Bam files and QC by running RNA-seq-Flow: `03A_Adzibolosu_fw_gg.sh` and `03B_EGA_8.25.25.sh`(Clauset_ABNexus/RNAseq_flow_out/)
  * Some EGA fastqs were of technical replicates. The bams corresponding to these files were therefore merged using `03BB_EGA_merge_bams.sh`
2. Count over the genes (all isoforms or the probe-focused exons) using `04A_Ad_gtf_featurecounts.sh` and `04B_EGA_gtf_featurecounts.sh` (Clauset_ABNexus/Counting/)
3. Clean up the counts: Cleans up the counts: links NanoString to isoforms, cleans up naming and gene names with `05_Get_Final_RNAseq_Counts.ipynb`

### Initial Comparison and calculation of L2FCs (Modeling)
1. `00_Compare_L2FC.ipynb`: Calculate log2FC, and do initial exploratory comparisons (e.g. PCA, range comparison)

### Prediction of NanoString vs RNAseq samples or PFS (Modeling)
Feature selection was done with two approaches: 
* Consistent: Only consider genes which were in the top X predictors across Y experiments (based on F-score).
* Bootstrapping: Only consider genes that were found in the top X predictors (based on F-score) across a fraction of at least X of bootstrapped samples (with weighted bootstrapping).

0. The functions used for predictive modeling, splitting up data, and graphing can be found at `helper_functions.py`
1. `01_NanoRNAseq_Cons.ipynb`: Modeling using Consistent feature selection, predicting NanoString vs RNA-seq, and PFS category.
2. `02_NanoRNAseq_Btsp.ipynb`: Modeling using Boostrapping feature selection, predicting NanoString vs RNA-seq, and PFS category.
3. `03_RNAseqFull_BtspCons.ipynb`: Feature selection using both consistent & bootstrapping approaches. Model prediction of PFS considering full RNA-seq data (both the normal RNA-seq dataset and using the more independent set)
