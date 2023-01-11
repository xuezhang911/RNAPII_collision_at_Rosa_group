#!/bin/bash -l
#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 72:00:00
#SBATCH -J RNAseq_alignment
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools star/2.7.9a

for i in $(find /crex/proj/snic2021-23-14/Xue/PRJDB10113/trimmomatic/  -name "*_trimmomatic.fq.gz");
do
a1=$(echo $i | sed -e 's/_trimmomatic.fq.gz//')
printf "this_is\t$a1\n"

STAR --genomeDir /crex/proj/snic2021-23-14/Xue/genome/arab_STAR_genome \

--readFilesIn $i \
--readFilesCommand zcat \ # we need to unzip the trimmomatic gz file \
--outFileNamePrefix /crex/proj/snic2021-23-14/Xue/PRJDB10113/RNA-seq_align_out/$a1 \
--limitBAMsortRAM 16000000000 \
--outSAMtype BAM SortedByCoordinate \ # convert sam file to bam
--outFilterType BySJout \
--outBAMsortingThreadN 5 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--quantMode TranscriptomeSAM GeneCounts
Done



