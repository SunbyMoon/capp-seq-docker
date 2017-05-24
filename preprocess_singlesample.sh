#!/bin/bash

#Arguments 
sample=$1 #CCD333-T-DNA
fastq=$2 #/data/s3/averafastq/patients/CCD333
out=$3 #/data/storage/capp-seq/patients/CCD333
repo=$4 #/home/anu/capp-seq-docker

#Variables
envlist=$out/$sample.env.list

# env.list
echo "sample="$sample >> $envlist
echo "fastq_dir="$fastq >> $envlist 
echo "aln_dir="$out"/alignment" >> $envlist 
echo "var_dir="$out"/vcf" >> $envlist 
echo "fastqc_dir="$out"/fastqc" >> $envlist
echo "bbduk_dir="$out"/bbduk" >> $envlist
echo "stats_dir="$out"/stats" >> $envlist
echo -e "\n" >> $envlist 

echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
echo "temp_dir=/tmp" >> $envlist
echo "regions=/home/panceq/data/PANCeq_C200X50bp_6.bed" >> $envlist #image
echo "regions=/data/database/avera/MedExome_hg19_capture_targets_noann.bed" >> $envlist
echo "annodb=/data/database/annovar" >> $envlist
echo -e "\n" >> $envlist 

echo "infastq1="$fastq"/"$sample"_1.fastq.gz" >> $envlist
echo "infastq2="$fastq"/"$sample"_2.fastq.gz" >> $envlist
echo "RGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$sample >> $envlist
echo "tmpsortbam=/tmp/"$sample"_sorted.bam" >> $envlist
echo "tmpdedupbam=/tmp/"$sample"_dedup.bam" >> $envlist
echo "tmprealnbam=/tmp/"$sample"_realigned.bam" >> $envlist
echo "tmpabra=/tmp/"$sample"_abra" >> $envlist
echo "tmpabralog=/tmp/"$sample"_abra.log" >> $envlist
echo "outbam="$out"/alignment/"$sample"_realigned_sorted.bam" >> $envlist
echo "fastq1="$out"/bbduk/"$sample"_1.fastq.gz" >> $envlist
echo "fastq2="$out"/bbduk/"$sample"_2.fastq.gz" >> $envlist
echo "bhist="$out"/stats/"$sample".bhist" >> $envlist
echo "qhist="$out"/stats/"$sample".qhist" >> $envlist
echo "aqhist="$out"/stats/"$sample".aqhist" >> $envlist
echo "lhist="$out"/stats/"$sample".lhist" >> $envlist
echo "gchist="$out"/stats/"$sample".gchist" >> $envlist
echo "qm_dir="$out"/qualimap" >> $envlist
echo -e "\n" >> $envlist 

echo "outvardict="$out"/vcf/"$sample"_vardict.vcf" >> $envlist
echo "outvardictgz="$out"/vcf/"$sample"_vardict.vcf.gz" >> $envlist
echo -e "\n" >> $envlist

mkdir -p $out/alignment
mkdir -p $out/vcf
mkdir -p $out/qualimap
mkdir -p $out/fastqc
mkdir -p $out/bbduk
mkdir -p $out/stats

sed -e 's|'/home/anu/capp-seq-docker/env.list'|'$envlist'|g' $repo/panceq_singlesample.sh > $out/$sample.panceq.sh
#bash $out/$sample.panceq.sh

