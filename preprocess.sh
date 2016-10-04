#!/bin/bash

#Arguments 
tsample=$1
gsample=$2
fastq=$3
out=$4
envlist=$5

sample=${tsample%%-*}
echo "tsample="$tsample >> $envlist
echo "gsample="$gsample >> $envlist 
echo "fastq_dir="$fastq >> $envlist 
echo "aln_dir="$out"/alignment" >> $envlist 
echo "var_dir="$out"/vcf" >> $envlist 
echo "cnv_dir="$out"/cnv" >> $envlist
echo "qm_dir="$out"/qualimap" >> $envlist
echo "fastqc_dir="$out"/fastqc" >> $envlist
echo -e "\n" >> $envlist 

echo "sample="${tsample%%-*} >> $envlist
echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
echo "temp_dir=/tmp" >> $envlist
echo "regions=/data/database/avera/CAPPSeq_Avera_capture_targets_6.bed" >> $envlist #image
echo "cont=/home/anu/capp-seq-docker/cont.R" >> $envlist #image
echo "contpanel=/home/anu/capp-seq-docker/contPanel.csv" >> $envlist #image
echo "cnvkit_filter=/home/anu/capp-seq-docker/cnvkit_filter.R" >> $envlist #image
echo -e "\n" >> $envlist 

echo "tfastq1="$fastq"/"$tsample"_1.fastq.gz" >> $envlist
echo "tfastq2="$fastq"/"$tsample"_2.fastq.gz" >> $envlist
echo "tRGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$tsample >> $envlist
echo "ttmpsortbam=/tmp/"$tsample"_sorted.bam" >> $envlist
echo "ttmpdedupbam=/tmp/"$tsample"_dedup.bam" >> $envlist
echo "ttmprealnbam=/tmp/"$tsample"_realigned.bam" >> $envlist
echo "ttmpabra=/tmp/"$tsample"_abra" >> $envlist
echo "ttmpabralog=/tmp/"$tsample"_abra.log" >> $envlist
echo "toutbam="$aln_dir/$tsample"_realigned_sorted.bam" >> $envlist
echo "tpileup="$aln_dir/$tsample".flt.pileup" >> $envlist #varscan
echo -e "\n" >> $envlist 

echo "gfastq1="$fastq"/"$gsample"_1.fastq.gz" >> $envlist
echo "gfastq2="$fastq"/"$gsample"_2.fastq.gz" >> $envlist
echo "gRGR=@RG\tID:1\tLB:SRR\tPL:ILLUMINA\tSM:"$gsample >> $envlist
echo "gtmpsortbam=/tmp/"$gsample"_sorted.bam" >> $envlist
echo "gtmpdedupbam=/tmp/"$gsample"_dedup.bam" >> $envlist
echo "gtmprealnbam=/tmp/"$gsample"_realigned.bam" >> $envlist
echo "gtmpabra=/tmp/"$gsample"_abra" >> $envlist
echo "gtmpabralog=/tmp/"$gsample"_abra.log" >> $envlist
echo "goutbam="$aln_dir/$gsample"_realigned_sorted.bam" >> $envlist
echo "gpileup="$aln_dir/$gsample".flt.pileup" >> $envlist #varscan
echo -e "\n" >> $envlist 

echo "outsom="$var_dir/$sample"_varscan_somatic" >> $envlist #varscan
echo "outvcfsnp="$var_dir/$sample"_varscan_somatic.snp.vcf" >> $envlist #varscan
echo "outvcfindel="$var_dir/$sample"_varscan_somatic.indel.vcf" >> $envlist #varscan
echo "outvcfsnphc="$var_dir/$sample"_varscan_somatic.snp.Somatic.hc.vcf" >> $envlist #varscan
echo "outvcfindelhc="$var_dir/$sample"_varscan_somatic.indel.Somatic.hc.vcf" >> $envlist #varscan
echo "outvcf="$var_dir/$sample"_varscan_somatic.vcf" >> $envlist #varscan
echo "outvcfgz="$var_dir/$sample"_varscan_somatic.vcf.gz" >> $envlist #varscan
echo "outvardict="$var_dir/$sample"_vardict_somatic.vcf.gz" >> $envlist
echo -e "\n" >> $envlist

echo "refcnvkit="$cnv_dir/$sample"_ref.cnn"
echo "cnvkitcns=/"$cnv_dir/$sample".cns"
echo "cnvkitout="$cnv_dir/$sample".cnvkit.out"

echo "contout="$var_dir/$sample".contamination.out" >> $envlist
echo "tqmout="$tsample"_qualimap"
echo "gqmout="$gsample"_qualimap"

mkdir -p $out/alignment
mkdir -p $out/vcf
mkdir -p $out/cnv
mkdir -p $out/qualimap
mkdir -p $out/fastqc

#sed -e 's|'/home/anu/capp-seq-docker/test.env.list'|'$envlist'|g' panceq.sh
#bash panceq.sh

