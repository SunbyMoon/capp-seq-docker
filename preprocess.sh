#!/bin/bash

#Arguments 
tsample=$1 #CCD333-T-DNA
gsample=$2 #CCD333-B-DNA
fastq=$3 #/data/s3/averafastq/patients/CCD333
out=$4 #/data/storage/capp-seq/patients/CCD333
repo=$5 #/home/anu/capp-seq-docker

#Variables
sample=${tsample%%-*}
id_tumor=${tsample#*-}
id_blood=${gsample#*-}
envlist=$out/$sample.env.list
qc_config=$out/$sample.yaml

# env.list
echo "tsample="$tsample >> $envlist
echo "gsample="$gsample >> $envlist 
echo "fastq_dir="$fastq >> $envlist 
echo "aln_dir="$out"/alignment" >> $envlist 
echo "var_dir="$out"/vcf" >> $envlist 
echo "cnv_dir="$out"/cnv" >> $envlist
echo "qm_dir="$out"/qualimap" >> $envlist
echo "fastqc_dir="$out"/fastqc" >> $envlist
echo "bbduk_dir="$out"/bbduk" >> $envlist
echo "stats_dir="$out"/stats" >> $envlist
echo "report_dir="$out"/report" >> $envlist
echo -e "\n" >> $envlist 

echo "sample="${tsample%%-*} >> $envlist
echo "id_tumor="$id_tumor >> $envlist
echo "id_blood="$id_blood >> $envlist 
echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
echo "temp_dir=/tmp" >> $envlist
echo "regions=/data/database/avera/CAPPSeq_Avera_capture_targets_6.bed" >> $envlist #image
echo "annodb=/data/database/annovar" >> $envlist
echo -e "\n" >> $envlist 

echo "tinfastq1="$fastq"/"$tsample"_1.fastq.gz" >> $envlist
echo "tinfastq2="$fastq"/"$tsample"_2.fastq.gz" >> $envlist
echo "tRGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$tsample >> $envlist
echo "ttmpsortbam=/tmp/"$tsample"_sorted.bam" >> $envlist
echo "ttmpdedupbam=/tmp/"$tsample"_dedup.bam" >> $envlist
echo "ttmprealnbam=/tmp/"$tsample"_realigned.bam" >> $envlist
echo "ttmpabra=/tmp/"$tsample"_abra" >> $envlist
echo "ttmpabralog=/tmp/"$tsample"_abra.log" >> $envlist
echo "toutbam="$out"/alignment/"$tsample"_realigned_sorted.bam" >> $envlist
echo "tfastq1="$out"/bbduk/"$tsample"_1.fastq.gz" >> $envlist
echo "tfastq2="$out"/bbduk/"$tsample"_2.fastq.gz" >> $envlist
echo "tbhist="$out"/stats/"$tsample".bhist" >> $envlist
echo "tqhist="$out"/stats/"$tsample".qhist" >> $envlist
echo "taqhist="$out"/stats/"$tsample".aqhist" >> $envlist
echo "tlhist="$out"/stats/"$tsample".lhist" >> $envlist
echo "tgchist="$out"/stats/"$tsample".gchist" >> $envlist
echo "tqm_dir="$out"/qualimap/tumor" >> $envlist
echo -e "\n" >> $envlist 

echo "ginfastq1="$fastq"/"$gsample"_1.fastq.gz" >> $envlist
echo "ginfastq2="$fastq"/"$gsample"_2.fastq.gz" >> $envlist
echo "gRGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$gsample >> $envlist
echo "gtmpsortbam=/tmp/"$gsample"_sorted.bam" >> $envlist
echo "gtmpdedupbam=/tmp/"$gsample"_dedup.bam" >> $envlist
echo "gtmprealnbam=/tmp/"$gsample"_realigned.bam" >> $envlist
echo "gtmpabra=/tmp/"$gsample"_abra" >> $envlist
echo "gtmpabralog=/tmp/"$gsample"_abra.log" >> $envlist
echo "goutbam="$out"/alignment/"$gsample"_realigned_sorted.bam" >> $envlist
echo "gfastq1="$out"/bbduk/"$gsample"_1.fastq.gz" >> $envlist
echo "gfastq2="$out"/bbduk/"$gsample"_2.fastq.gz" >> $envlist
echo "gbhist="$out"/stats/"$gsample".bhist" >> $envlist
echo "gqhist="$out"/stats/"$gsample".qhist" >> $envlist
echo "gaqhist="$out"/stats/"$gsample".aqhist" >> $envlist
echo "glhist="$out"/stats/"$gsample".lhist" >> $envlist
echo "ggchist="$out"/stats/"$gsample".gchist" >> $envlist
echo "gqm_dir="$out"/qualimap/normal" >> $envlist
echo -e "\n" >> $envlist 

echo "outvardict="$out"/vcf/"$sample"_vardict_somatic.vcf" >> $envlist
echo "outvardictgz="$out"/vcf/"$sample"_vardict_somatic.vcf.gz" >> $envlist
echo "annovarout="$out"/vcf/"$sample"_vardict_somatic.annovar" >> $envlist
echo "annovcf="$out"/vcf/"$sample"_vardict_somatic.annovar.hg19_multianno.vcf" >> $envlist
echo "filtvcf="$out"/vcf/"$sample"_somatic_postfilter.vcf" >> $envlist
echo -e "\n" >> $envlist

echo "refcnvkit="$out"/cnv/"$sample"_ref.cnn" >> $envlist
echo "cnvkitcns="$out"/cnv/"$tsample"_realigned_sorted.cns" >> $envlist
echo "cnvkitout="$out"/cnv/"$sample".cnvkit.out" >> $envlist
echo -e "\n" >> $envlist

echo "contout="$out"/report/"$sample".contamination.out" >> $envlist
echo "qc_config="$qc_config >> $envlist

mkdir -p $out/alignment
mkdir -p $out/vcf
mkdir -p $out/cnv
mkdir -p $out/qualimap/tumor
mkdir -p $out/qualimap/normal
mkdir -p $out/fastqc
mkdir -p $out/bbduk
mkdir -p $out/stats
mkdir -p $out/report

sed -e 's|'/home/anu/capp-seq-docker/env.list'|'$envlist'|g' $repo/panceq.sh > $out/$sample.panceq.sh
sed -e "1s/.*/$sample:/" $repo/sample_test.yaml > $out/$sample.yaml
#bash $out/$sample.panceq.sh

