# Step1 : build up conda environment in terminal
conda  info --envis
Conda activate xue
# download package blast in conda environment
Conda install -c bioconda blast
# install perl
conda install perl-digest-md5
# download wget function
# set up brew tool https://code2care.org/howto/install-homebrew-brew-on-m1-mac
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
Brew install wget
# download arabidopsis fasta file
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 
gzip -d Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 
# make a blast database
Makeblastdb -in Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -dbtype nucl -out TAIR10 -parse_seqids
## modify the path of blast 
# no idea why, but creat a .ncbirc file at home
Touch .ncbirc
BLASTDB=/User/xue/Database/blast
BLASTDB_NUCL_DATA_LOADER=/User/xue/Database/blast/nt
WINDOW_MASKER_PATH=/User/xue/Database/blast/windowmasker
# use update_blastdb.pl to update and download nt file ) nt.00.tar.gz  nt.00.tar.gz.md5
nohup time update_blastdb.pl nt nr > log &
# creat a local fasta file
Nano q.fa
# for example the complementary sequence corresponding to one probe sequence from SOFL1 TTGAAAAAGAAGAGCATAGG
# blast
Blastn -task blastn -db TAIR10 -query q.fa -outfmt 7 -out q.txt
Head -n 15 q.txt
Blastn -task blastn -db TAIR10 -query q.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qseq sseq evalue bitscore" -out q.txt
wc q.txt
less q.txt
# filter q.txt with the overlap length greater than 16bp
awk -F "\t" '$4-$5-$6 >16' q.txt>nq.txt
# dowload arabidopsis gtf file 
wget ftp://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.54.gtf.gz
mv Arabidopsis_thaliana.TAIR10.54.gtf ath.gtf
# extract the chromosome sequence
cat  ath.gtf  | perl -alne '{next unless $F[2] eq "gene" 
# or optionally : sed 's/"/\t/g' ath.gtf |awk -v type="gene" 'BEGIN{OFS=FS="\t"}{if($3=="gene") {print $1,$4-1,$5,$10,".",$7}}' > ath.geneid.bed
#  install bedtools
 conda install -c bioconda bedtools
#  prepared bed file 
Mv nq.txt nq.bed
cat nq.bed |awk -F "\t" '{if($9<$10) {print $2,$9-1,$10}}'>1.bed
cat nq.bed |awk -F "\t" '{if($9>$10) {print $2,$10-1,$9}}'>2.bed
cat 1.bed >>2.bed
# check if the file is tab delimited 
 cat -t 2.bed
#make it as a tab delimited file 
perl -p -i -e 's/ /\t/g' 2.bed

# generate a file with gene name 
bedtools intersect -a 2.bed -b ath.geneid.bed -wa -wb >a.bed 
# if we add chr to the column 1, we can do sed 's/^/chr&/g' file >file1

