
#!/bin/bash

#Arguments
tsample=$1 #CCD333-T-DNA
gsample=$2 #CCD333-B-DNA
fastq=$3 #/data/s3/averafastq/patients/CCD333
out=$4 #/data/storage/capp-seq/patients/CCD333
repo=$5 #/home/anu/capp-seq-docker
bed=$6 #/data/database/avera/PANCeq_C200X50bp_6.bed
cnvbed=$7 #/data/database/avera/PANCeq_CNV_capture_targets_6.bed
genome=$8 #/data/database/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa
anno=$9 #/home/anu/repos/annovar
annodb=$10 #/data/database/annovar
oncodb=$11 #/data/database/avera/oncotator_v1_ds_April052016
temp=$12 #/tmp


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

echo "tumor: "$tsample
echo "normal: "$gsample
echo "fastq_dir: "$fastq
echo "out_dir: "$out
echo "repo_dir: "$repo
#echo "bed: "$bed
#echo "cnvbed: "$cnvbed
#echo "genome: "$genome
#echo "anno: "$anno
#echo "annodb: "$annodb
#echo "oncodb: "$oncodb
#echo "temp: "$temp


#if [ -z ${bed+x} ]; then echo "bed is unset"; else echo "bed is set to '$bed'"; fi
#if [ -z ${cnvbed+x} ]; then echo "cnvbed is unset"; else echo "cnvbed is set to '$cnvbed'"; fi
#if [ -z ${genome+x} ]; then echo "genome is unset"; else echo "genome is set to '$genome'"; fi
#if [ -z ${anno+x} ]; then echo "anno is unset"; else echo "anno is set to '$anno'"; fi
#if [ -z ${annodb+x} ]; then echo "annodb is unset"; else echo "annodb is set to '$annodb'"; fi
#if [ -z ${oncodb+x} ]; then echo "oncodb is unset"; else echo "oncodb is set to '$oncodb'"; fi
#if [ -z ${temp+x} ]; then echo "temp is unset"; else echo "temp is set to '$temp'"; fi
