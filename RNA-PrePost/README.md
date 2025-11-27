## Paired HGSOC Pre & Post NACT Sample Analysis

### Organizing the Metadata and NanoString Counts

Used to get NonEGA, and NanoString formatted Counts

1. Get all metadata into one file and format NanoString counts. Nanostring/analysis/`Get_metadata_and_format_counts.ipynb`
2. EGA data had several sources of metadata that needed to be linked together. This was done using `Harmonizing_Metadata.ipynb` --> should be EGA 
3. Combine all metadata with `Combine_EGA_Rest_Metadata.ipynb`

### Get RNA-seq Counts
Javella et al. only provided already normalized counts so we could not map and count for ourselves.
1. Get the Bam files and QC by running RNA-seq-Flow: `Adzibolosu_fw_gg.sh` and `EGA_8.25.25.sh`(Clauset_ABNexus/RNAseq_flow_out/)
  * Some EGA fastqs were of technical replicates. The bams corresponding to these files were therefore merged using `EGA_merge_bams.sh`
2. Count over the genes (all isoforms or the probe-focused exons) using `Ad_gtf_featurecounts.sh` and `EGA_gtf_featurecounts.sh` (Clauset_ABNexus/Counting/)
3. Clean up the counts: Cleans up the counts: links NanoString to isoforms, cleans up naming and gene names with `Get_Final_RNAseq_Counts.ipynb`

### Initial Comparison and calculation of L2FCs (Exploratory)
1. `Compare_L2FC_11.24.25.ipynb`: Calculate log2FC, and do initial exploratory comparisons

### Prediction of NanoString vs RNAseq samples or PFS
Feature selection was done with two approaches: 
* Consistent: Only consider genes which were in the top X predictors across Y experiments (based on F-score).
* Bootstrapping: Only consider genes that were found in the top X predictors (based on F-score) across a fraction of at least X of bootstrapped samples (with weighted bootstrapping).

1. `Cons_Iter_Feature_11.3.25.ipynb`: Feature selection using consistent approaches for model prediction of NanoString and PFS.
2. `BootStrap_Feature_11.3.25.ipynb`: Feature selection using bootstrapping approaches for model prediction of NanoString and PFS.
3. `ModelPred_FullRnaIsof_11.7.25.ipynb`: Feature selection using both consistent & bootstrapping approaches. Model prediction of PFS considering full RNA-seq data (both the normal RNA-seq dataset and using the more independent set)