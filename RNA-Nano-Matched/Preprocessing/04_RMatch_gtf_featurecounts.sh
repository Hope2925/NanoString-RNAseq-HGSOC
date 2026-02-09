#!/bin/bash
#SBATCH --job-name=Ad_Counting
#SBATCH --output=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/RMatch_featurecounts_%j.out
#SBATCH --error=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/RMatch_featurecounts_%j.out
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=5G
#SBATCH --partition short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hoto7260@colorado.edu

##########################
# EDIT THE FOLLOWING
##########################

#### Loading Modules/Environments

module load bedtools
module load samtools/1.8
module load subread/1.6.2

BAM_DIR=/scratch/Shares/clauset/Clauset_ABNexus/RNAseq_flow_out/RNA_Matched_4.18/mapped/bams
GTF=/scratch/Shares/clauset/Clauset_ABNexus/Counting/gtf_counts/data/Refseq_GRCh38.p14.gtf
EXON_GTF=/scratch/Shares/clauset/Clauset_ABNexus/Counting/gtf_counts/data/PanCancer_IO_360_Probe_Exon_hg38.gtf
COUNT_DIR=/scratch/Shares/clauset/Clauset_ABNexus/Counting/gtf_counts/out
 
########################################## 
# Count reads over gene coordinates     ##
##########################################
cd ${bams}
echo "Submitted counts with Rsubread for original genes by strand......"
# Use this many threads (-T)
# Count multi-overlapping read
# Count reads in a strand specific manner (-s 2)
# Count by the gene feature (-t 'gene')
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -t "exon" \
    -a ${GTF} \
    -F 'GTF' \
    -o ${COUNT_DIR}/RNAMatch4.18_str_gtf_genes.txt \
    ${BAM_DIR}/*.sorted.bam

featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -t "exon" \
    -a ${GTF} \
    -F 'GTF' \
    -g transcript_id \
    -o ${COUNT_DIR}/RNAMatch4.18_str_gtf_transcripts.txt \
    ${BAM_DIR}/*.sorted.bam


########################################## 
# Count reads over Exon coordinates     ##
##########################################
cd ${bams}
echo "Submitted counts with Rsubread for original genes by strand......"
# Use this many threads (-T)
# Count multi-overlapping read
# Count reads in a strand specific manner (-s 2)
# Count by the gene feature (-t 'gene')
featureCounts \
    -T 32 \
    -O \
    -s 1 \
    -t "exon" \
    -a ${EXON_GTF} \
    -F 'GTF' \
    -o ${COUNT_DIR}/RNAMatch4.18_str_gtf_matchedexons.txt \
    ${BAM_DIR}/*.sorted.bam

