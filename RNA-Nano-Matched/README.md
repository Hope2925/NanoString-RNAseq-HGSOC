## Matched Data

### Prep the data
1. Ran RNA-seq-flow with `RNA_Matched_3.19.25.sh` on the fastq files (Clauset_ABNexus/RNAseq_flow_out). This got us the bams.
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
2. `Evaluate_Harmonization_Post.ipynb`