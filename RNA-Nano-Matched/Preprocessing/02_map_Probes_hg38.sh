#!/bin/bash
#SBATCH --job-name=HISAT2 # Job name
#SBATCH -p short
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hope.townsend@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=8     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=30gb # Memory limit
#SBATCH --time=6:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/Hisat2map.%j.out # Standard output
#SBATCH --error=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/Hisat2map.%j.out # Standard error log

# Get the files
HISAT2_INX=/scratch/Shares/dowell/genomes/hg38/HISAT2/genome
SEQ_FILE=/scratch/Shares/clauset/Clauset_ABNexus//Counting/Probe_prep/data/PanCancer_IO_360_Probe_Sequences.fastq
OUT_SAM=/scratch/Shares/clauset/Clauset_ABNexus/Counting/Probe_prep/out/PanCancer_IO_360_Probe_hg38_hisat.sam
OUT_STATS=/scratch/Shares/clauset/Clauset_ABNexus/Counting/Probe_prep/out/PanCancer_IO_360_Probe_hg38.hisat2_mapstats.txt  
# Load modules
module load hisat2/2.1.0

# Max and min penalties for soft clipping and mismatch: 1,0 & 3,1

hisat2  -p 8 \
                --very-sensitive \
                --rna-strandness F \
                --pen-noncansplice 14 \
                --mp 1,0 \
                --sp 3,1 \
                -x ${HISAT2_INX} \
                -U ${SEQ_FILE} \
                --new-summary \
                > ${OUT_SAM} \
                2> ${OUT_STATS}   


