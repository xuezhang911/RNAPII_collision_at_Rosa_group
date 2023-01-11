# we use fastqc to perform quality control of fasta sequence
#!/bin/bash -l

#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 36:00:00
#SBATCH -J  fastqc
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

module load bioifanfo-tools FastQC/0.11.9
mkdir fastqc
fastqc -o  
/crex/proj/snic2021-23-14/Xue/tissue_RNA_seq/fastqc   /crex/proj/snic2021-23-14/Xue/tissue_RNA_seq/raw/*.fastq.gz

