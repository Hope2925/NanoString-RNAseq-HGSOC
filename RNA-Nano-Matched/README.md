## Matched Data

Scripts and notebooks are found under notebooks/

### Prep the data
1. Ran [RNA-seq-flow](https://github.com/Dowell-Lab/RNAseq-Flow) with `RNA_Matched_3.19.25.sh` on the fastq files (available in GEO X) (Clauset_ABNexus/RNAseq_flow_out). This got us the bams and QC results. 
### Match Probes to Isoforms/Exons
* Clauset_ABNexus/Probe_pre/bin/
1. `Get_Probes_fasta.ipynb` creates fastq sequences from the probe regions that will be used for HISAT2
  * Input: PanCancer_IO_360_Gene_List.txt
2. `map_Probes_hg38.sh` maps the probes to hg38 using HISAT2
3. `Get_Probe_Exon_GTF.ipynb`
  * Overlaps probe coordinates according to HISAT2 to exons, gets the coordinates of the exons and the converts to a GTF for easy counting 
  * Output: PanCancer_IO_360_Probe_Exon_hg38.gtf

### Counting
1. `Rmatch_gtf_featurecounts.sh` (Clauset_ABNexus/Counting/)
  * Gets the counts for the exons, isoforms, and genes
2. `Get_Final_Counts_10.1.25.ipynb`
  * Cleans up the counts: links NanoString to isoforms, cleans up naming and gene names

### Final Comparisons
1. `Evaluate_Harmonization_10.20.25.ipynb`
  * Evaluate correlation of NanoString and RNAseq expression (across 6 harmonization approaches) both across patients (per-gene coefficients) and within patients (per-patient coefficients)
  * Show what trends (e.g. expression levels or exon number) might explain why some genes have poor correlation, regardless of harmonization approach
    * Estimate limits of detection based on expression levels and range of expression for which NanoString and RNA-seq can be fairly compared
  * Show what trends (e.g. exon number ) might explain why some genes have better harmonization when considering different harmonization approaches (e.g. Isoform- vs Exon-based counts, or length-based vs other normalization approaches)
2. `Evaluate_Harmonization_Post.ipynb`
  * Evaluate impact on harmonization from removing cases below limits of detection.