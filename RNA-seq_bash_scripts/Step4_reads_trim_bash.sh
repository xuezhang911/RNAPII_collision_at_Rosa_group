
#We use trimmomatic remove adapters. The name of bash file is trimmomatic_bash.sh
#!/bin/bash -l
sample=$(basename ${1/_sortmerna.fq.gz})
echo "sample name is $sample"
module load bioinfo-tools trimmomatic/0.36 FastQC/0.11.9  MultiQC/1.11
trimmomatic SE $1 /crex/proj/snic2021-23-14/Xue/PRJDB10113/trimmomatic/${sample}_trimmomatic.fq.gz ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50
fastqc -o /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/fastqc-trimmomatic --noextract -t 1 /crex/proj/snic2021-23-14/Xue/PRJDB10113/trimmomatic/${sample}_trimmomatic.fq.gz
multiqc -f -o  /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/fastqc-trimmomatic/ .

