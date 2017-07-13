#!/bin/bash

#Arguments 
tsample=$1 #CCD333-T-DNA
gsample=$2 #CCD333-B-DNA
fastq=$3 #/data/s3/averafastq/patients/CCD333
out=$4 #/data/storage/capp-seq/patients/CCD333
repo=$5 #/home/anu/capp-seq-docker

bed=/data/database/avera/MedExome_hg19_capture_targets_noann_6col.bed
#bed=/data/database/avera/PANCeq_C200X50bp_6.bed
cnvbed=/data/database/avera/MedExome_hg19_capture_targets_noann_6col.bed
#cnvbed=/data/database/avera/PANCeq_CNV_capture_targets_6.bed
genome=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
anno=/home/anu/repos/annovar
annodb=/data/database/annovar
oncodb=/data/database/avera/oncotator_v1_ds_April052016
temp=/tmp

#bed=$6 #/data/database/avera/PANCeq_C200X50bp_6.bed
#cnvbed=$7 #/data/database/avera/PANCeq_CNV_capture_targets_6.bed
#genome=$8 #/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
#anno=$9 #/home/anu/capp-seq-docker
#annodb=$10 #/data/database/annovar
#oncodb=$11 #/data/database/avera/oncotator_v1_ds_April052016
#temp=$12 #/tmp

# parse arguments
usage () {
echo -e "\033[1mUsage:\033[0m bash preprocess.sh -t TUMOR_NAME -n NORMAL_NAME -f FASTQ_DIR -o OUT_DIR -r REPOSITORY [-b BED] [-c CNV_BED] [-g GENOME] [-a ANNOVAR_DIR] [-d ANNOVAR_DB] [-o ONCOTATOR_DB] [-t TEMP_DIR] [-h]"
echo -e "\n"
echo -e "\033[1mRequired options:\033[0m"
echo "  -t TUMOR_NAME, --tumor TUMOR_NAME		Name of the tumor sample. Ex: CCD333-T-DNA"
echo "  -n NORMAL_NAME, --normal NORMAL_NAME	Name of the normal sample. Ex: CCD333-B-DNA"
echo "  -f FASTQ_DIR, --fastq FASTQ_DIR		Full path of the directory where FASTQ files are located"
echo "						(tumor and normal FASTQ files are expected to be in the same location)"
echo "						Ex: /data/s3/averafastq/patients/CCD333"
echo "  -o OUT_DIR, --out OUT_DIR			Full path of the output directory"
echo "						Subfolders will be created within this directory"
echo "						Ex: /data/storage/capp-seq/patients/CCD333"
echo "  -r REPOSITORY, --repo REPOSITORY 		Full path of the PANCeq repository. Ex: /home/anu/capp-seq/docker"
echo -e "\n"
echo -e "\033[1mOptional arguments:\033[0m"
echo "  -b BED, --bed BED				Full path for the BED file to be used for alignment and variant calling"
echo "						Default: /data/database/avera/PANCeq_C200X50bp_6.bed"
echo "  -c CNV_BED, --cnvbed CNV_BED 			Full path of the BED file to be used for CNV calling"
echo "						Default: /data/database/avera/PANCeq_CNV_capture_targets_6.bed"
echo "  -g GENOME, --genome GENOME 			Full path of the genome FASTA file"
echo "						Default: /data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
echo "  -a ANNOVAR_DIR, --anno ANNOVAR_DIR 		Full path of the directory where ANNOVAR scripts are stored"
echo "						Default: /home/anu/repos/annovar"
echo "  -d ANNOVAR_DB, --annodb ANNOVAR_DB		Full path of the ANNOVAR database"
echo "						Default: /data/database/annovar"
echo "  -n ONCOTATOR_DB, oncodb ONCOTATOR_DB 		Full path of the ONCOTATOR database"
echo "						Default: /data/database/avera/oncotator_v1_ds_April052016"
echo "  -p TEMP_DIR, --temp TEMP_DIR 			Full path of the temporary directory to be used"
echo "						Default: /tmp"
echo "  -h, --help 					Program usage is displayed"
}

#options=`getopt -o t:n:f:o:r:b::c::g::a::d::n::p::h --l tumor:,normal:,fastq:,out:,repo:,bed::,cnvbed::,genome::,anno::,annodb::,oncodb::,temp::,help -- "$@"`
options=`getopt -o t:n:f:o:r:h --l tumor:,normal:,fastq:,out:,repo:,help -- "$@"`
eval set -- "$options"
#echo $options

while true; do
	case "${1}" in
		-t|--tumor)
			tsample=${2} ; shift 2 ;;
		-n|--normal)
			gsample=${2} ; shift 2 ;;
		-f|--fastq)
			fastq=${2} ; shift 2 ;;
		-o|--out)
			out=${2} ; shift 2 ;;
		-r|--repo)
			repo=${2} ; shift 2 ;;
#		-b|--bed)
#			case "${2}" in
#				"") bed=/data/database/avera/PANCeq_C200X50bp_6.bed ; shift 2 ;;
#				*) bed=${2} ; shift 2 ;;
#			esac ;;
#                -c|--cnvbed)
#                        case "${2}" in
#                                "") cnvbed=/data/database/avera/PANCeq_CNV_capture_targets_6.bed ; shift 2 ;;
#                                *) cnvbed=${2} ; shift 2 ;;
#			esac ;;
#                -g|--genome)
#                        case "${2}" in
#                                "") genome=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa ; shift 2 ;;
#                                *) genome=${2} ; shift 2 ;;
#			esac ;;
#                -a|--anno)
#                        case "${2}" in
#                                "") anno=/home/anu/repos/annovar ; shift 2 ;;
#                                *) anno=${2} ; shift 2 ;;
#			esac ;;
#                -d|--annodb)
#                        case "${2}" in
#                                "") annodb=/data/database/annovar ; shift 2 ;;
#                                *) annodb=${2} ; shift 2 ;;
#			esac ;;
#                -n|--oncodb)
#                        case "${2}" in
#                                "") oncodb=/data/database/avera/oncotator_v1_ds_April052016 ; shift 2 ;;
#                                *) oncodb=${2} ; shift 2 ;;
#			esac ;;
#                -p|--temp)
#                        case "${2}" in
#                                "") temp=/tmp ; shift 2 ;;
#                                *) temp=${2} ; shift 2 ;;
#			esac ;;
		-h|--help)
			usage
			exit;;
                --) shift ; break ;;
		*) echo "Internal error!" ; exit 1 ;;
	esac
done

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
echo "repo="$repo >> $envlist
echo -e "\n" >> $envlist 

echo "sample="${tsample%%-*} >> $envlist
echo "id_tumor="$id_tumor >> $envlist
echo "id_blood="$id_blood >> $envlist 
echo "cpu="$(grep -c "processor" /proc/cpuinfo) >> $envlist
echo -e "\n" >> $envlist

echo "regions="$bed >> $envlist
echo "cnvregions="$cnvbed >> $envlist
echo "bwa_index="$genome >> $envlist
echo "anno="$anno >> $envlist
echo "annodb="$annodb >> $envlist
echo "oncodb="$oncodb >> $envlist
echo "temp_dir="$temp >> $envlist
echo -e "\n" >> $envlist

#echo "bwa_index=/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa" >> $envlist
#echo "temp_dir=/tmp" >> $envlist
#echo "regions=/data/database/avera/PANCeq_C200X50bp_6.bed" >> $envlist #image
#echo "regions=/data/database/avera/MedExome_hg19_capture_targets_noann.bed" >> $envlist
#echo "annodb=/data/database/annovar" >> $envlist

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

echo "oncoout="$out"/vcf/"$sample"_vardict_somatic.oncotator.vcf" >> $envlist
echo "oncooutgz="$out"/vcf/"$sample"_vardict_somatic.oncotator.vcf.gz" >> $envlist
echo -e "\n" >> $envlist

echo "refcnvkit="$out"/cnv/"$sample"_ref.cnn" >> $envlist
echo "cnvkitcns="$out"/cnv/"$tsample"_realigned_sorted.cns" >> $envlist
echo "cnvkitout="$out"/cnv/"$sample".cnvkit.out" >> $envlist
echo -e "\n" >> $envlist

echo "tmbout="$out"/report/"$sample".tmb.out" >> $envlist
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

