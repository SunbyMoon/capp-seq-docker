#!/bin/bash

#Arguments 
sample=$1 #CCD333-T-DNA
#gsample=$2 #CCD333-B-DNA
fastq=$2 #/data/s3/averafastq/patients/CCD333
out=$3 #/data/storage/capp-seq/patients/CCD333
repo=$4 #/home/anu/capp-seq-docker

#Variables
#sample=${tsample%%-*}
#id_tumor=${tsample#*-}
#id_blood=${gsample#*-}
envlist=$out/$sample.env.list
#qc_config=$out/$sample.yaml

# env.list
echo "sample="$sample >> $envlist
#echo "gsample="$gsample >> $envlist 
echo "fastq_dir="$fastq >> $envlist 
echo "aln_dir="$out"/alignment" >> $envlist 
echo "var_dir="$out"/vcf" >> $envlist 
#echo "cnv_dir="$out"/cnv" >> $envlist
#echo "qm_dir="$out"/qualimap" >> $envlist
echo "fastqc_dir="$out"/fastqc" >> $envlist
echo "bbduk_dir="$out"/bbduk" >> $envlist
echo "stats_dir="$out"/stats" >> $envlist
#echo "report_dir="$out"/report" >> $envlist
echo -e "\n" >> $envlist 

#echo "sample="${tsample%%-*} >> $envlist
#echo "id_tumor="$id_tumor >> $envlist
#echo "id_blood="$id_blood >> $envlist 
echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
echo "temp_dir=/tmp" >> $envlist
echo "regions=/data/database/avera/PANCeq_C200X50bp_6.bed" >> $envlist #image
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

#echo "ginfastq1="$fastq"/"$gsample"_1.fastq.gz" >> $envlist
#echo "ginfastq2="$fastq"/"$gsample"_2.fastq.gz" >> $envlist
#echo "gRGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$gsample >> $envlist
#echo "gtmpsortbam=/tmp/"$gsample"_sorted.bam" >> $envlist
#echo "gtmpdedupbam=/tmp/"$gsample"_dedup.bam" >> $envlist
#echo "gtmprealnbam=/tmp/"$gsample"_realigned.bam" >> $envlist
#echo "gtmpabra=/tmp/"$gsample"_abra" >> $envlist
#echo "gtmpabralog=/tmp/"$gsample"_abra.log" >> $envlist
#echo "goutbam="$out"/alignment/"$gsample"_realigned_sorted.bam" >> $envlist
#echo "gfastq1="$out"/bbduk/"$gsample"_1.fastq.gz" >> $envlist
#echo "gfastq2="$out"/bbduk/"$gsample"_2.fastq.gz" >> $envlist
#echo "gbhist="$out"/stats/"$gsample".bhist" >> $envlist
#echo "gqhist="$out"/stats/"$gsample".qhist" >> $envlist
#echo "gaqhist="$out"/stats/"$gsample".aqhist" >> $envlist
#echo "glhist="$out"/stats/"$gsample".lhist" >> $envlist
#echo "ggchist="$out"/stats/"$gsample".gchist" >> $envlist
#echo "gqm_dir="$out"/qualimap/normal" >> $envlist
#echo -e "\n" >> $envlist 

echo "outvardict="$out"/vcf/"$sample"_vardict.vcf" >> $envlist
echo "outvardictgz="$out"/vcf/"$sample"_vardict.vcf.gz" >> $envlist
#echo "annovarout="$out"/vcf/"$sample"_vardict_somatic.annovar" >> $envlist
#echo "annovcf="$out"/vcf/"$sample"_vardict_somatic.annovar.hg19_multianno.vcf" >> $envlist
#echo "filtvcf="$out"/vcf/"$sample"_somatic_postfilter.vcf" >> $envlist
echo -e "\n" >> $envlist

#echo "refcnvkit="$out"/cnv/"$sample"_ref.cnn" >> $envlist
#echo "cnvkitcns="$out"/cnv/"$tsample"_realigned_sorted.cns" >> $envlist
#echo "cnvkitout="$out"/cnv/"$sample".cnvkit.out" >> $envlist
#echo -e "\n" >> $envlist

#echo "contout="$out"/report/"$sample".contamination.out" >> $envlist
#echo "qc_config="$qc_config >> $envlist

mkdir -p $out/alignment
mkdir -p $out/vcf
#mkdir -p $out/cnv
mkdir -p $out/qualimap
#mkdir -p $out/qualimap/normal
mkdir -p $out/fastqc
mkdir -p $out/bbduk
mkdir -p $out/stats
#mkdir -p $out/report

sed -e 's|'/home/anu/capp-seq-docker/env.list'|'$envlist'|g' $repo/panceq_singlesample.sh > $out/$sample.panceq.sh
#sed -e "1s/.*/$sample:/" $repo/sample_test.yaml > $out/$sample.yaml
#bash $out/$sample.panceq.sh

