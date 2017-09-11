#!/bin/bash

#===============================================================================
#
#          FILE:  neoantigen_pipeline_vcf.sh
#
# 
#   DESCRIPTION:  This is an pipeline calling mutation from vcf file
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   ChiZhou 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  5/15/2017 15:09 CST
#      REVISION:  ---
#===============================================================================

#######--Argument--######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: ChiZhou"
	echo "	 This is an pipeline calling mutation from vcf file"
	echo "   Usage:"
	echo "	 bash neoantigen_pipeline_vcf.sh [option] <-i input_vcf> <-o out_dir> <-P prefix> <-t HLA_type> <-r reference> <-c cpu_num> <-b bindingaffinity_cutoff> <-E expression_file> <-f fpkm_cutoff>"
	echo ""
	echo "   Required argument:"
	echo "  	-i input vcf file. mutation in vcf file format "
	echo "  	-t 4-digit HLA genotype,eg:HLA-A02:01,if more than one,seperate by comma"
	echo "   Optional arguments:"
	echo "  	-P prefix for output file. default: vcf file name"
	echo "  	-o output directory. default: ./"
	echo "  	-c CPU thread. --default: 4"
	echo "  	-C site reads coverge cutoff. default: 200"
	echo " 		-a binding affinity fold change cutoff .default: 1"
	echo "  	-b binding affinity cutoff .default: 500"
	echo "  	-E expression profile of tumor sample"
	echo "  	-f gene expression fpkm value cutoff, this option should not set when -E is empty"
	echo ""
	echo -e "\033[31m  !!! NOTE: \033[0m"
	echo -e "	The program can deal with pair-end tumor wgs data. If not, please modify the program."
	echo "	Program needed: bwa edwstates picards samtools bedtools Varscan VEP HLAminer netMHC netchop Sequenza PyClone "
	echo "	Make sure all the program is contained by environment viariable PATH."
	echo "	The pipeline also need index files builded already, if not the pipeline will build it."
	echo ""
	echo "    (-^V^-) Enjoy yourself~~~"
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo ""
}


#OUTPATH=./
#PREFIX=0
#CPU=4
#Binding_Aff_Fc_Cutoff=1
#Binding_Aff_Cutoff=500
#Exp_file=""
#Fpkm_Cutoff=1
if [ $# -lt 2 ];then
	help_info
	exit 1
fi

TEMP=`getopt -o i:t:o:P:c:a:b:E:f: \
--long input_vcf:,HLA_type:,outpath:,prefix:,CPU:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff: \
	-n 'neoantigen_pipeline_vcf.sh' -- "$@"`

if [ $? != 0 ];then
	echo "Terminating..." >&2 ; exit 1 ; fi


eval set -- "$TEMP"
while true
do
	case "$1" in
		-i | --input_vcf) vcf_file=$2;shift 2;;
		-t | --HLA_type) HLA_type_str=$2;shift 2;;
		-o | --outpath) OUTPATH=$2;shift 2;;
		-P | --prefix) PREFIX=$2;shift 2;;
		-c | --CPU) CPU=$2;shift 2;;
		-b | --binding_aff_cutoff) Binding_Aff_Cutoff=$2;shift 2;;
		-a | --binding_aff_fc_cutoff) Binding_Aff_Fc_Cutoff=$2;shift 2;;
		-E | --expression_file) Exp_file=$2;shift 2;;
		-f | --fpkm_cutoff) Fpkm_Cutoff=$2;shift 2;;
		--) shift;break;;
		*) echo "Internal error!";exit 1;;
	esac

done
###########check Argument setting###########
if [ ! -f "${vcf_file}" ];then
	echo -e "\033[40;31;1mERROR: please use -i to specify the right WGS fastq file\033[0m"
	exit 1
fi

if [ ! -n "${HLA_type_str}" ];then
	echo -e "\033[40;31;1mERROR: please use -t to specify the HLA type of sample\033[0m"
	exit 1
fi
echo ${vcf_file}
echo ${HLA_type_str}

if [ ! -n "$OUTPATH" ];then
	OUTPATH=./${PREFIX}
	mkdir ./${PREFIX}
	echo -e "\033[40;33;1mWARNING: output path not found, use current fold\033[0m"
elif [ ! -d "$OUTPATH/${PREFIX}" ];then
	echo -e "\033[40;33;1mWARNING: output path not found, create one\033[0m"	
	mkdir -p ${OUTPATH}/${PREFIX}
fi


if [ ! -n "$PREFIX" ];then
	echo -e "\033[40;33;1mWARNING: no prefix name, use sample name as prefix\033[0m"
	PREFIX=${N_fastq_1%.f*q*}
	PREFIX=${PREFIX##*/}
fi

if [ ! -n "$Binding_Aff_Cutoff" ];then
	Binding_Aff_Cutoff=500
#elif [ $Binding_Aff_Cutoff -lt 0 ];then
#	echo "the Binding_Aff_Cutoff should be greater than 0!"
#	exit 1
fi

#if [ "$Exp_file" = "None" ];then
#	Exp_file="no_exp"
#fi

if [ ! -d "$Exp_file" ];then
	Exp_file="no_exp"
fi


if [ ! -n "$Binding_Aff_Fc_Cutoff" ];then
	Binding_Aff_Fc_Cutoff=1
#elif [ $Binding_Aff_Fc_Cutoff -lt 0 ];then
#	echo "the Binding_Aff_Fc_Cutoff should be greater than 0!"
#	exit 1
fi
if [ ! -n "$Fpkm_Cutoff" ];then
	Fpkm_Cutoff=2
#elif [ $Fpkm_Cutoff -lt 0 ];then
#	echo "the Fpkm_Cutoff should be greater than 0!"
#	exit 1
fi

######result fold prepparation#######
######make result directories######

DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Preparation......\t"$DATE"\033[0m"
cd ${OUTPATH}/${PREFIX}

if [ ! -d liftover ];then
	mkdir liftover
fi

if [ ! -d somatic_mutation ];then
	mkdir somatic_mutation
fi

if [ ! -d netmhc ];then
	mkdir netmhc
fi

if [ ! -d netctl ];then
	mkdir netctl
fi

if [ ! -d annotation ];then
	mkdir annotation
fi


####liftover hg19 coordinates to hg38###
#CrossMap.py vcf /home/zhouchi/database/Annotation/crossmap/hg19ToHg38.over.chain.gz ${vcf_file} /home/zhouchi/database/Annotation/Fasta/human.fasta liftover/${PREFIX}_liftover38.vcf
###split vcf file into snv and indel#####
#SNP
#vcftools --vcf liftover/${PREFIX}_liftover38.vcf --remove-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}_SNVs_only
#INDEL
#vcftools --vcf liftover/${PREFIX}_liftover38.vcf --keep-only-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}_INDELs_only

####not liftover####
#SNP
#vcftools --gzvcf ${vcf_file} --remove-indels --recode --recode-INFO-all --stdout | gzip -c > somatic_mutation/${PREFIX}_SNVs_only.vcf.gz
#INDEL
#vcftools --gzvcf ${vcf_file} --keep-only-indels --recode --recode-INFO-all --stdout | gzip -c > somatic_mutation/${PREFIX}_INDELs_only.vcf.gz

netmhcpan(){
	local input_fasta=$1
	local hla_str=$2
	local netmhc_out=$3
	local out_dir=$4
	local split_num=$5
	if [ -d ${out_dir}/tmp ];then
		rm -rf ${out_dir}/tmp
		mkdir -p ${out_dir}/tmp
	else
		mkdir -p ${out_dir}/tmp
	fi
	if [ -f ${out_dir}/${netmhc_out} ];then
		rm ${out_dir}/${netmhc_out}
	fi
	split -l ${split_num} ${input_fasta} ${out_dir}/tmp/
	filelist=`ls ${out_dir}/tmp/`
	arr1=(${filelist})
	echo ${arr1[@]}
	OLD_IFS="$IFS" 
	IFS=","
	arr2=(${hla_str})
	IFS="$OLD_IFS" 
	for s in ${arr2[@]}
	do
	{
		echo $s
		for file_l in ${arr1[@]}
		do
		{
			echo ${file_l}
			netMHCpan -a $s -f ${out_dir}/tmp/${file_l} -l 8,9,10,11 > ${out_dir}/tmp/${s}_${file_l}_tmp_netmhc.txt
		} &
		done
		wait
		
	}
	done
	for file_l in ${arr1[@]}
	do
	{
		rm ${out_dir}/tmp/${file_l}
	}
	done
	filelist1=`ls ${out_dir}/tmp/`
	for file_r in $filelist1
	do
	{
		cat ${out_dir}/tmp/${file_r} >> ${out_dir}/${netmhc_out}
		rm ${out_dir}/tmp/${file_r}	
	}
	done
	rm -rf 	${out_dir}/tmp
}


####SNV derived neoantigens####
echo -e "\033[40;32mVEP Annotation......\t"$DATE"\033[0m"
VEP_CACHE=/home/zhouchi/database/Annotation/vep_data_89_37
#vep.pl -i ${vcf_file} --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_snv_vep_ann.txt --force_overwrite
#python ${iTuNES_BIN_PATH}/snv2fasta.py -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -o netmhc/ -s ${PREFIX}
echo -e "\033[40;32mnetMHCpan......\t"$DATE"\033[0m"
set -x
#netmhcpan netmhc/${PREFIX}_snv.fasta ${HLA_type_str} ${PREFIX}_snv_netmhc.txt netmhc 20000
set +x
#######filtering######
set -x
#python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_snv_netmhc.txt -g netmhc/${PREFIX}_snv.fasta -o netmhc -s ${PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${HLA_type_str}
set +x
#######netchop########
set -x
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_final_neo_candidate.txt -o netctl -s ${PREFIX}
#python ${iTuNES_BIN_PATH}/neoantigens_annotation.py -i netctl/${SNV_PREFIX}_netctl_concact.txt -o annotation -s ${PREFIX}






####INDEL derived neoantigens###
set -x
#vep.pl -i somatic_mutation/${PREFIX}_INDELs_only.vcf.gz --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_indel_vep_ann.txt --force_overwrite
#gunzip somatic_mutation/${PREFIX}_INDELs_only.vcf.gz
#python ${iTuNES_BIN_PATH}/deletion2fasta.py -i somatic_mutation/${PREFIX}_indel_vep_ann.txt -c somatic_mutation/${PREFIX}_INDELs_only.vcf -o netmhc -s ${PREFIX}
#python ${iTuNES_BIN_PATH}/insertion2fasta.py -i somatic_mutation/${PREFIX}_indel_vep_ann.txt -c somatic_mutation/${PREFIX}_INDELs_only.recode.vcf -o netmhc -s ${PREFIX}
set +x
########netMHC###########
#cat netmhc/${PREFIX}_del.fasta > netmhc/${PREFIX}_indel.fasta
#cat netmhc/${PREFIX}_del.fasta > netmhc/${PREFIX}_indel.fasta
#netmhcpan netmhc/${PREFIX}_indel.fasta ${hla_str} ${PREFIX}_snv_netmhc.txt netmhc 100
#######filtering######
set -x
#python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_indel_netmhc.txt -g netmhc/${PREFIX}_indel.fasta -o netmhc -s ${PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
set +x
#######netchop########
set -x
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${INDEL_PREFIX}_final_neo_candidate.txt -o netctl -s ${PREFIX}
set +x





















