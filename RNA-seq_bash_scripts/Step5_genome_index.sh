# gtf file and genome file are required
# we index Arabidopsis thalina genome with STAR
#!/bin/bash -l
#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 72:00:00
#SBATCH -J genomeindex
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools star/2.7.9a
STAR --runMode genomeGenerate \
--genomeDir /crex/proj/snic2021-23-14/Xue/genome \
--genomeFastaFiles /crex/proj/snic2021-23-14/Xue/genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--sjdbGTFfile /crex/proj/snic2021-23-14/Xue/genome/ara.gtf \
--sjdbOverhang 149

