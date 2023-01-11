# create a loop for running sortmerna software
#!/bin/bash -l

#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 12:00:00
#SBATCH -J sortmerna
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL


dir /crex/proj/snic2021-23-14/Xue/PRJDB10113/rawreads/
for i in $(seq 54 58);
do bash /crex/proj/snic2021-23-14/Xue/scripts/sortmerna_bash.sh /crex/proj/snic2021-23-14/Xue/PRJDB10113/rawreads/DRR2353$i.fastq.gz
done

