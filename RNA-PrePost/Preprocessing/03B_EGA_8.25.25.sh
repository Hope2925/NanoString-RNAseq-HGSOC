#!/bin/bash
#SBATCH --job-name=RNAnextflow # Job name
#SBATCH -p long
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hope.townsend@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=8gb # Memory limit
#SBATCH --time=200:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/EGA_8.27.25_RNAflow.%j.out # Standard output
#SBATCH --error=/scratch/Shares/clauset/Clauset_ABNexus/e_and_o/EGA_8.27.25_RNAflow.%j.err # Standard error log

#activate virtual environment
source /scratch/Shares/clauset/Clauset_ABNexus/venv/bin/activate
echo "Activated environment"
#load modules 
module load sra/2.8.0
module load bbmap/38.05
module load fastqc/0.11.8
module load hisat2/2.1.0
module load samtools/1.8
module load preseq
module load igvtools/2.3.75
module load mpich/3.2.1
module load bedtools/2.28.0
module load openmpi/1.6.4
module load gcc/7.1.0
module load python/3.6.3/rseqc/3.0.0


#Get the Nextflow paths (for script and Nextflow Executive)
SRC=/Users/hoto7260/Flows/RNAseq-Flow
NF_EXE='/scratch/Shares/dowell/dbnascent/pipeline_assets/nextflow'


${NF_EXE} ${SRC}/main.nf -profile slurm --workdir '/scratch/Shares/clauset/Clauset_ABNexus/RNAseq_flow_out/tmp/' --genome_id 'hg38' --outdir '/scratch/Shares/clauset/Clauset_ABNexus/RNAseq_flow_out/EGA_8.27/' --email hope.townsend@colorado.edu --fastqs '/scratch/Shares/clauset/Clauset_ABNexus/EGA/fastqs/*.fastq.gz' --reverseStranded

## They originally note the 1 as Forward and _2 as Reverse



