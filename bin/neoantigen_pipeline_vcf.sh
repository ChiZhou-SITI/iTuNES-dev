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
	OUTPATH=./
	echo -e "\033[40;33;1mWARNING: output path not found, use current fold\033[0m"
elif [ ! -d "$OUTPATH" ];then
	echo -e "\033[40;33;1mWARNING: output path not found, create one\033[0m"	
	mkdir -p ${OUTPATH}
fi
echo $OUTPATH

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
cd ${OUTPATH}

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

if [ ! -d liftover/${PREFIX} ];then
	mkdir liftover/${PREFIX}
fi

if [ ! -d somatic_mutation/${PREFIX} ];then
	mkdir somatic_mutation/${PREFIX}
fi

if [ ! -d netmhc/${PREFIX} ];then
	mkdir netmhc/${PREFIX}
fi
if [ ! -d netctl/${PREFIX} ];then
	mkdir netctl/${PREFIX}
fi
if [ ! -d annotation/${PREFIX} ];then
	mkdir annotation/${PREFIX}
fi

####liftover hg19 coordinates to hg38###
#CrossMap.py vcf /home/zhouchi/database/Annotation/crossmap/hg19ToHg38.over.chain.gz ${vcf_file} /home/zhouchi/database/Annotation/Fasta/human.fasta liftover/${PREFIX}/${PREFIX}_liftover38.vcf


###split vcf file into snv and indel#####
#SNP
#vcftools --vcf liftover/${PREFIX}/${PREFIX}_liftover38.vcf --remove-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}/${PREFIX}_SNVs_only
#INDEL
#vcftools --vcf liftover/${PREFIX}/${PREFIX}_liftover38.vcf --keep-only-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}/${PREFIX}_INDELs_only



####SNV derived neoantigens####
echo -e "\033[40;32mVEP Annotation......\t"$DATE"\033[0m"
VEP_CACHE=/data/PUBLIC/VEP_DATA
#vep.pl -i somatic_mutation/${PREFIX}/${PREFIX}_SNVs_only.recode.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}/${PREFIX}_snv_vep_ann.txt --force_overwrite
SNV_PREFIX=${PREFIX}_snv
#python ${iTuNES_BIN_PATH}/snv2fasta.py -i somatic_mutation/${PREFIX}/${PREFIX}_snv_vep_ann.txt -o netmhc/${PREFIX} -s ${PREFIX}
########netMHC###########
#cd netmhc/${PREFIX}
echo -e "\033[40;32mnetMHCpan......\t"$DATE"\033[0m"
#if [ -d tmp ];then
#	rm -rf tmp/
#else
#	mkdir tmp
#fi
#if [ -f ${PREFIX}_snv_netmhc.txt ];then
#	rm ${PREFIX}_snv_netmhc.txt
#fi
#cd tmp
#split -l 1000 ../${PREFIX}_snv.fasta tmp
#filelist=`ls .`
#for file in $filelist
#do
#{
#	netMHCpan -a ${HLA_type_str} -f ${file} -l 9 > ${file}_tmp_netmhc.txt
#	rm ${file}
#} &
#done
#wait
#filelist1=`ls .`
#for file_r in $filelist1
#do
#{
#	cat $file_r >> ../${PREFIX}_snv_netmhc.txt
#	rm ${file_r}	
#}
#done
#cd ..
#rm -rf tmp
#cd ..
#cd ..
#######filtering######
set -x
#python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}/${PREFIX}_snv_netmhc.txt -g netmhc/${PREFIX}/${PREFIX}_snv.fasta -o netmhc/${PREFIX} -s ${SNV_PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
set +x
#######netchop########
set -x
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}/${SNV_PREFIX}_final_neo_candidate.txt -o netctl/${PREFIX} -s ${SNV_PREFIX}
#python ${iTuNES_BIN_PATH}/neoantigens_annotation.py -i netctl/${PREFIX}/${SNV_PREFIX}_netctl_concact.txt -o annotation/${PREFIX} -s ${SNV_PREFIX}





####INDEL derived neoantigens###
#vep.pl -i somatic_mutation/${PREFIX}/${PREFIX}_INDELs_only.recode.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}/${PREFIX}_indel_vep_ann.txt --force_overwrite
INDEL_PREFIX=${PREFIX}_indel
#python ${iTuNES_BIN_PATH}/deletion2fasta.py -i somatic_mutation/${PREFIX}/${PREFIX}_indel_vep_ann.txt -c somatic_mutation/${PREFIX}/${PREFIX}_INDELs_only.recode.vcf -o netmhc/${PREFIX} -s ${INDEL_PREFIX}
########netMHC###########
set -x
#cd netmhc/${PREFIX}
#echo -e "\033[40;32mnetMHCpan......\t"$DATE"\033[0m"
#if [ -d tmp ];then
#	rm -rf tmp/
#else
#	mkdir tmp
#fi
#if [ -f ${PREFIX}_indel_netmhc.txt ];then
#	rm ${PREFIX}_indel_netmhc.txt
#fi
#cd tmp
#split -l 1000 ../${PREFIX}_indel.fasta tmp
#filelist=`ls .`
#for file in $filelist
#do
#{
#	netMHCpan -a $HLA_type_str -f ${file} -l 9 > ${file}_tmp_netmhc.txt
#	rm ${file}
#} &
#done
#wait
#filelist1=`ls .`
#for file_r in $filelist1
#do
#{
#	cat $file_r >> ../${PREFIX}_indel_netmhc.txt
#	rm ${file_r}	
#}
#done
#cd ..
#rm -rf tmp
#cd ..
#cd ..
set +x
#######filtering######
set -x
#python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}/${PREFIX}_indel_netmhc.txt -g netmhc/${PREFIX}/${PREFIX}_indel.fasta -o netmhc/${PREFIX} -s ${INDEL_PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
set +x
#######netchop########
set -x
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}/${INDEL_PREFIX}_final_neo_candidate.txt -o netctl/${PREFIX} -s ${INDEL_PREFIX}
#python ${iTuNES_BIN_PATH}/neoantigens_annotation.py -i netclt/${PREFIX}/${INDEL_PREFIX}_netctl_concact.txt -o annotation/${PREFIX} -s ${INDEL_PREFIX}






















