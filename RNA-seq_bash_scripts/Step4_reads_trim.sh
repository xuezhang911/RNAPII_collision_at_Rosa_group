#create a loop for running trimmomatic software
#!/bin/bash -l
#SBATCH -A snic2021-22-12
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 30:00:00
#SBATCH -J trimmomatic
#SBATCH --mail-user xue.zhang@slu.se
#SBATCH --mail-type=ALL
dir /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/
for fwd in $(find /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/  -name "*_sortmerna.fq.gz");
do bash /crex/proj/snic2021-23-14/Xue/scripts/trimmomatic_bash.sh $fwd
done

