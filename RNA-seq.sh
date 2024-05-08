#!/bin/bash

mkdir rnaseq && cd rnaseq
mkdir rawdata
mkdir ref

############################ Download reference and annotation #######################################
cd ref
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
gunzip gencode.v36.annotation.gtf.gz
cd ..



##################################### fastqc ########################################
mkdir outputs
mkdir outputs/fastqc
fastqc ./rawdata/*fastq.gz -o ./outputs/fastqc -t 2



##################################### multiqc ###########################################
mkdir outputs/multiqc
multiqc ./outputs/fastqc -o ./outputs/multiqc


################################# Filter low quality sequences ###############################
mkdir ./outputs/fastqc_trimmed

# 提取所有样本名称
samples=$(ls rawdata/ | grep '_R1.fastq' | awk -F '_R1.fastq' '{print $1}')
# echo查看样本名称
echo $samples

for sample in $samples
do
    trim_galore --fastqc -j 4 --paired --basename $sample -o ./outputs/fastqc_trimmed "rawdata/"$sample"_R1.fastq.gz" "rawdata/"$sample"_R2.fastq.gz"
done 



###################################### Build index ###################################
mkdir outputs/hisat2_index
hisat2_extract_splice_sites.py ref/gencode.v36.annotation.gtf > outputs/hisat2_index/gencode.v36.annotation.ss
hisat2_extract_exons.py ref/gencode.v36.annotation.gtf > outputs/hisat2_index/gencode.v36.annotation.exon
hisat2-build -p 4 --ss outputs/hisat2_index/gencode.v36.annotation.ss --exon outputs/hisat2_index/gencode.v36.annotation.exon ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna outputs/hisat2_index/gencode.v36.annotation_tran





######################################## Alignment ##################################
mkdir ./outputs/hisat2_alignment

for sample in $samples
do
    hisat2 -p 4 --dta -x outputs/hisat2_index/gencode.v36.annotation_tran -1 "outputs/fastqc_trimmed/"$sample"_R1_val_1.fq.gz" -2 "outputs/fastqc_trimmed/"$sample"_R2_val_2.fq.gz" -S "outputs/hisat2_alignment/"$sample".sam"
done



#################################### sam to bam ##################################
mkdir ./outputs/samtools_bam

for sample in $samples
do
    samtools view -bS "outputs/hisat2_alignment/"$sample".sam" > "outputs/samtools_bam/"$sample".bam"
    samtools sort "outputs/samtools_bam/"$sample".bam" -o "outputs/samtools_bam/"$sample".sort.bam"
done


################################# featureCounts #######################################
mkdir ./outputs/featureCounts

for sample in $samples
do
    featureCounts -p -T 4 -t exon -g gene_id -a ref/gencode.v36.annotation.gtf -o "outputs/featureCounts/"$sample"_count.tsv" "outputs/samtools_bam/"$sample".sort.bam"
done