#!/bin/bash -l
#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -J getreads
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2021-23-14/Xue/raw
module load bioinfo-tools sratools/2.9.6-1
for i in $(seq 28 32);
do fastq-dump --gzip --split-3  DRR2353$i
# if it's pair-end
gzip DRR2353$i"_1.fastq" DRR2353$i"_2.fastq" 
#if it's single_end: gzip DRR2353$i".fastq"
done
