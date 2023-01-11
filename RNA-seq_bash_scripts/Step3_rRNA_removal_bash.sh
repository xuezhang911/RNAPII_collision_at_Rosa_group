#We use sortmeRNA remove rRNA. The name of bash file is sortmerna_bash.sh
#!/bin/bash -l
set -eu
if [ -d /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/kvdb ]; then
 rm -rf /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/kvdb/*
Fi

sample=$(basename ${1/.fastq.gz})
echo "sample name is $sample"
module load bioinfo-tools SortMeRNA/4.3.3
sortmerna --ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5.8s-database-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/rfam-5s-database-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-16s-id95.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-arc-23s-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-16s-id90.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-bac-23s-id98.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-18s-id95.fasta \
--ref /crex/proj/snic2021-23-14/Xue/rRNArefgenome/silva-euk-28s-id98.fasta \
--reads $1 \
--workdir /crex/proj/snic2021-23-14/Xue/sortmerna1 --fastx --other 
mv /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/out/aligned.log /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/out/${sample}_aligned.log
mv /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/out/other.fq  /crex/proj/snic2021-23-14/Xue/PRJDB10113/sortmerna/out/${sample}_sortmerna.fq

