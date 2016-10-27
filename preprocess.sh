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
echo "bbduk_dir="$out"/bbduk" >> $envlist
echo "stats_dir="$out"/stats" >> $envlist
echo -e "\n" >> $envlist 

echo "sample="${tsample%%-*} >> $envlist
echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
echo "temp_dir=/tmp" >> $envlist
echo "regions=/data/database/avera/CAPPSeq_Avera_capture_targets_6.bed" >> $envlist #image
echo "cont=/home/anu/capp-seq-docker/cont.R" >> $envlist #image
echo "contpanel=/home/anu/capp-seq-docker/contPanel.csv" >> $envlist #image
echo "cnvkit_filter=/home/anu/capp-seq-docker/cnvkit_filter.R" >> $envlist #image
echo "cnvkit_regions=/home/anu/capp-seq-docker/PANCeq_CNV_capture_targets_6.bed" >> $envlist #image
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
echo "tpileup="$out"/alignment/"$tsample".flt.pileup" >> $envlist #varscan
echo "tfastq1="$out"/bbduk/"$tsample"_1.fastq.gz" >> $envlist
echo "tfastq2="$out"/bbduk/"$tsample"_2.fastq.gz" >> $envlist
echo "tbhist="$out"/stats/"$tsample".bhist" >> $envlist
echo "tqhist="$out"/stats/"$tsample".qhist" >> $envlist
echo "taqhist="$out"/stats/"$tsample".aqhist" >> $envlist
echo "tlhist="$out"/stats/"$tsample".lhist" >> $envlist
echo "tgchist="$out"/stats/"$tsample".gchist" >> $envlist
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
echo "gpileup="$out"/alignment/"$gsample".flt.pileup" >> $envlist #varscan
echo "gfastq1="$out"/bbduk/"$gsample"_1.fastq.gz" >> $envlist
echo "gfastq2="$out"/bbduk/"$gsample"_2.fastq.gz" >> $envlist
echo "gbhist="$out"/stats/"$gsample".bhist" >> $envlist
echo "gqhist="$out"/stats/"$gsample".qhist" >> $envlist
echo "gaqhist="$out"/stats/"$gsample".aqhist" >> $envlist
echo "glhist="$out"/stats/"$gsample".lhist" >> $envlist
echo "ggchist="$out"/stats/"$gsample".gchist" >> $envlist
echo -e "\n" >> $envlist 

echo "outsom="$out"/vcf/"$sample"_varscan_somatic" >> $envlist #varscan
echo "outvcfsnp="$out"/vcf/"$sample"_varscan_somatic.snp.vcf" >> $envlist #varscan
echo "outvcfindel="$out"/vcf/"$sample"_varscan_somatic.indel.vcf" >> $envlist #varscan
echo "outvcfsnphc="$out"/vcf/"$sample"_varscan_somatic.snp.Somatic.hc.vcf" >> $envlist #varscan
echo "outvcfindelhc="$out"/vcf/"$sample"_varscan_somatic.indel.Somatic.hc.vcf" >> $envlist #varscan
echo "outvcf="$out"/vcf/"$sample"_varscan_somatic.vcf" >> $envlist #varscan
echo "outvcfgz="$out"/vcf/"$sample"_varscan_somatic.vcf.gz" >> $envlist #varscan
echo "outvardict="$out"/vcf/"$sample"_vardict_somatic.vcf.gz" >> $envlist
echo -e "\n" >> $envlist

echo "refcnvkit="$out"/cnv/"$sample"_ref.cnn" >> $envlist
echo "cnvkitcns="$out"/cnv/"$sample".cns" >> $envlist
echo "cnvkitout="$out"/cnv/"$sample".cnvkit.out" >> $envlist

echo "contout="$out/$sample".contamination.out" >> $envlist
echo "tqmout="$tsample"_qualimap" >> $envlist
echo "gqmout="$gsample"_qualimap" >> $envlist

mkdir -p $out/alignment
mkdir -p $out/vcf
mkdir -p $out/cnv
mkdir -p $out/qm
mkdir -p $out/fastqc
mkdir -p $out/bbduk
mkdir -p $out/stats

#sed -e 's|'/home/anu/capp-seq-docker/test.env.list'|'$envlist'|g' panceq.sh
#bash panceq.sh

