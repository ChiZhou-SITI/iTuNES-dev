#!/bin/bash

#===============================================================================
#
#          FILE:  neoantigen_pipeline_DNAnoMatched_se.sh
#
# 
#   DESCRIPTION:  This is an pipeline calling mutation from tumor WGSS raw data
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   ChiZhou 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  3/18/2017 22:013 CST
#      REVISION:  ---
#===============================================================================

#######--Argument--######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: ChiZhou"
	echo "	 This is an pipeline calling mutation from  single-end WGSS raw data"
	echo "   Usage:"
	echo "	 bash neoantigen_pipeline_DNAnoMatched_se.sh [option] <-p tumor.fq(fastq)> <-o out_dir> <-g genome> <-P prefix> <-I bwa_index> <-r reference> <-c cpu_num> <-C coverage> <-b bindingaffinity_cutoff> <-E expression_file> <-f fpkm_cutoff>"
	echo ""
	echo "   Required argument:"
	echo "  	-p tumor fastq1. first fastq file from tumor sample "
	echo "  	-g genome species.eg :human "
	
	echo "   Optional arguments:"
	echo "  	-I  bwa index. default: /home/zhouchi/database/Annotation/Index/bwa/genome"
	echo "  	-P prefix for output file. default: sample name"
	echo "  	-o output directory. default: ./"
	echo "  	-r genome reference file. default: /home/zhouchi/database/Annotation/Fasta/genome.fa"
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
#COVERAGE=200
#Binding_Aff_Fc_Cutoff=1
#Binding_Aff_Cutoff=500
#Exp_file=""
#Fpkm_Cutoff=1
if [ $# -lt 2 ];then
	help_info
	exit 1
fi

TEMP=`getopt -o p:o:g:P:I:r:c:C:a:b:E:f: \
--long T_fq:,outpath:,genome:,prefix:,bwa_index:,reference:,CPU:Coverage:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff: \
	-n 'neoantigen_pipeline_DNAMatched_se.sh' -- "$@"`

if [ $? != 0 ];then
	echo "Terminating..." >&2 ; exit 1 ; fi


eval set -- "$TEMP"
while true
do
	case "$1" in
		-p | --T_fq) T_fastq=$2;shift 2;;
		-o | --outpath) OUTPATH=$2;shift 2;;
		-g | --genome) GENOME=$2;shift 2;;
		-P | --prefix) PREFIX=$2;shift 2;;
		-I | --bwa_index) BWA_INDEX=$2;shift 2;;
		-r | --reference) REFERENCE=$2;shift 2;;
		-c | --CPU) CPU=$2;shift 2;;
		-C | --Coverage) COVERAGE=$2;shift 2;;
		-b | --binding_aff_cutoff) Binding_Aff_Cutoff=$2;shift 2;;
		-a | --binding_aff_fc_cutoff) Binding_Aff_Fc_Cutoff=$2;shift 2;;
		-E | --expression_file) Exp_file=$2;shift 2;;
		-f | --fpkm_cutoff) Fpkm_Cutoff=$2;shift 2;;
		--) shift;break;;
		*) echo "Internal error!";exit 1;;
	esac

done
###########check Argument setting###########
echo $T_fastq
if [ ! -f "$T_fastq" ];then
	echo -e "\033[40;31;1mERROR: please use -p to specify the right WGS fastq file\033[0m"
	exit 1
fi
if [ ! -n "$GENOME" ];then
	echo -e "\033[40;31;1mERROR: please use -g to specify genome used, eg. mm10\033[0m"
	exit 1
fi

if [ ! -n "$OUTPATH" ];then
	OUTPATH=./
	echo -e "\033[40;33;1mWARNING: output path not found, use current fold\033[0m"
elif [ ! -d "$OUTPATH" ];then
	echo -e "\033[40;33;1mWARNING: output path not found, create one\033[0m"	
	mkdir -p ${OUTPATH}
fi


if [ ! -n "$PREFIX" ];then
	echo -e "\033[40;33;1mWARNING: no prefix name, use sample name as prefix\033[0m"
	PREFIX=${N_fastq_1%.f*q*}
	PREFIX=${PREFIX##*/}
fi
#if [ ! -n "$BWA_INDEX" ];then
#	BWA_INDEX=/home/zhouchi/database/Annotation/Index/bwa/${GENOME}
#fi

#if [ ! -n "$REFERENCE" ];then
#	REFERENCE=/home/zhouchi/database/Annotation/Fasta/${GENOME}.fasta
#fi

if [ "$BWA_INDEX" = "None" ];then
	BWA_INDEX=/home/zhouchi/database/Annotation/Index/bwa/${GENOME}
fi

if [ "$REFERENCE" = "None" ];then
	REFERENCE=/home/zhouchi/database/Annotation/Fasta/${GENOME}.fasta
fi

if [ ! -n "$COVERAGE" ];then
	COVERAGE=100
#elif [ $COVERAGE -lt 0 ];then
#	echo "the COVERAGE should be greater than 0!"
#	exit 1
fi

if [ ! -n "$Binding_Aff_Cutoff" ];then
	Binding_Aff_Cutoff=500
#elif [ $Binding_Aff_Cutoff -lt 0 ];then
#	echo "the Binding_Aff_Cutoff should be greater than 0!"
#	exit 1
fi

if [ "$Exp_file" = "None" ];then
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


######result fold preparation#######
######make result directories######

DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Preparation......\t"$DATE"\033[0m"
cd ${OUTPATH}
if [ ! -d log_file ];then
	mkdir log_file
fi

if [ ! -d alignments ];then
	mkdir alignments
fi

if [ ! -d figures ];then
	mkdir figures
fi


if [ ! -d summary ];then
	mkdir summary
fi


if [ ! -d bamstat ];then
	mkdir bamstat
fi

if [ ! -d somatic_mutation ];then
	mkdir somatic_mutation
fi

if [ ! -d copynumber_profile ];then
	mkdir copynumber_profile
fi
if [ ! -d pyclone ];then
	mkdir pyclone
fi
if [ ! -d HLAtyping ];then
	mkdir HLAtyping
fi
if [ ! -d netmhc ];then
	mkdir netmhc
fi
if [ ! -d netctl ];then
	mkdir netctl
fi
echo ${GENOME} 
echo ${BWA_INDEX}
###check index
if [ -f ${BWA_INDEX}/${GENOME}".amb" -a -f ${BWA_INDEX}/${GENOME}".ann" -a -f ${BWA_INDEX}/${GENOME}".bwt" -a -f ${BWA_INDEX}/${GENOME}".pac" -a -f ${BWA_INDEX}/${GENOME}".sa" ];then
	echo -e "\033[40;35;1mIndex: "${BWA_INDEX}"\033[0m"
else
	echo "no index file, make bwa index"
	if [ -f ${GENOME_FA} ];then
		DATE=`date --date="-24 hour"`
		echo -e "\033[40;32mmake index......\t"$DATE"\033[0m"
		set -x
		bwa index -p ${BWA_INDEX} -a bwtsw ${REFERENCE}
		set +x
	else
		echo -e "\033[40;31;1mERROR: no index file and no genome.fa to build it\033[0m"
		exit 1
	fi
	echo -e "\033[40;35;1mIndex: "${BWA_INDEX}"\033[0m"
fi


#######check other files########





####################
###Align with bwa###
####################
###align
DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Align......---\t"$DATE"\033[0m"
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mbwa alignment......\t"$DATE"\033[0m"


######somatic mutation calling#########
Picard=/home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar
GATK="/home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar"
BuildBamIndex="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard BuildBamIndex"
MarkDuplicates="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard MarkDuplicates"
AddOrReplaceReadGroups="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard AddOrReplaceReadGroups"
RealignerTargetCreator="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T RealignerTargetCreator -nt 10 -dt NONE"
IndelRealigner="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T IndelRealigner"
BaseRecalibrator="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T BaseRecalibrator -nct 10"
AnalyzeCovariates="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T AnalyzeCovariates -dt NONE"
HaplotypeCaller="java -Djava.io.tmpdir=/home/zhouchi/tmp -Xmx8g -jar $GATK -T HaplotypeCaller -nct 10 -dt NONE"
SortSam="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard SortSam"
PrintReads="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T PrintReads -nct 10 -dt NONE"
VariantRecalibrator="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T VariantRecalibrator"
ApplyRecalibration="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T ApplyRecalibration"
VariantFiltration="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T VariantFiltration"
echo -e "\033[40;32mcalling somatic mutaiton using GATK......\t"$DATE"\033[0m"


BWA_sefq2bam(){
    local input_fq1=$1
    local output_bam=$2
    local output_dup=$3
    sam=`${iTuNES_BIN_PATH}/newtf`.sam;	#mkfifo $sam;
    bwa mem -t ${CPU} ${BWA_INDEX}/${GENOME} $input_fq1 > $sam                      # Make the output of bwa as bam
    sorted=`${iTuNES_BIN_PATH}/newtf`.bam;	#mkfifo $sorted;
    $SortSam I=$sam O=$sorted SO=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true                # MarkDuplicates cannot accept piped in sam file
    add_output_bam=`${iTuNES_BIN_PATH}/newtf`.bam; #mkfifo $add_output_bam;
    $AddOrReplaceReadGroups I=$sorted O=$add_output_bam SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20  # Add your own annotation
    $MarkDuplicates I=$add_output_bam O=$output_bam METRICS_FILE=$output_dup REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
    $BuildBamIndex I=$output_bam O=$output_bam.bai

    rm -f  $sam $sorted $add_output_bam
}

realignment(){
    local input_bam=$1
    local known_indels1=$2
    local known_indels2=$3
    local output=$4
    local intervals=`${iTuNES_BIN_PATH}/newtf`.intervals; #mkfifo $intervals;
    $RealignerTargetCreator -R $REFERENCE -I $input_bam -known $known_indels1 -known $known_indels2 -filterRNC -o $intervals 
    $IndelRealigner         -R $REFERENCE -I $input_bam -known $known_indels1 -known $known_indels2 -targetIntervals $intervals -filterRNC -o $output      # The parallel -nt -nct are not supported
    rm $intervals 
}

baserecalibration(){
    local input_bam=$1
    local known_indels1=$2
    local known_indels2=$3
    local known_indels3=$4
    local output=$5
    #local recal_plots=$6
    local before_table=`${iTuNES_BIN_PATH}/newtf`
    #local after_table=`${iTuNES_BIN_PATH}/newtf`
    $BaseRecalibrator  -R $REFERENCE -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3 -o $before_table
    $PrintReads        -R $REFERENCE -I $input_bam -BQSR $before_table -o $output
    #$BaseRecalibrator  -R $REFERENCE -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3  -BQSR $before_table -o $after_table
    #$AnalyzeCovariates -R $REFERENCE -before $before_table -after $after_table -plots $recal_plots -l DEBUG        
    rm $before_table
    # note that Rscript will be executed. many packages are needed. add -l DEBUG to check which library are needed 
    # if the execution failed. The following library are required:
    # reshape, gplots, ggplot2, bitops, caTools, colorspace, gtools, gData/Storage1, RColorBrewer, gsalib
}



callvariants_dna(){
    local input_dna_bam=$1
    local output=$2
    $HaplotypeCaller -R $REFERENCE -I $input_dna_bam -o $output -stand_call_conf 30 --dbsnp $dbsnp138
}

dbsnp138='/home/zhouchi/database/Annotation/hg38_vcf/dbsnp_138.hg38.vcf'
hapmap='/home/zhouchi/database/Annotation/hg38_vcf/hapmap_3.3.hg38.vcf'
omni='/home/zhouchi/database/Annotation/hg38_vcf/1000G_phase1.snps.high_confidence.hg38.vcf'
OneKG='/home/zhouchi/database/Annotation/hg38_vcf/1000G_omni2.5.hg38.vcf'
mills='/home/zhouchi/database/Annotation/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf'
set -x
#BWA_sefq2bam	$T_fastq	alignments/${PREFIX}_tumor_bwa.bam	alignments/${PREFIX}_tumor_bwa_dedup.dup
#realignment             alignments/${PREFIX}_tumor_bwa.bam        $OneKG     $mills  alignments/${PREFIX}_tumor_bwa_realign.bam
#baserecalibration       alignments/${PREFIX}_tumor_bwa_realign.bam      $dbsnp138      $OneKG     $mills   alignments/${PREFIX}_tumor_bwa_realign_recal.bam
#callvariants_dna            alignments/${PREFIX}_tumor_bwa_realign_recal.bam	somatic_mutation/${PREFIX}_raw.vcf
set +x
###split vcf file into snv and indel#####
#SNP
#vcftools --vcf somatic_mutation/${PREFIX}_raw.vcf --remove-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}/${PREFIX}_SNVs_only
#INDEL
#vcftools --vcf somatic_mutation/${PREFIX}_raw.vcf --keep-only-indels --recode --recode-INFO-all --out somatic_mutation/${PREFIX}/${PREFIX}_INDELs_only


##########neoantigen identification#####
#######HLA TYPING#######
:<<!
echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
#python $HLAtyping -i ${T_fastq_1} ${T_fastq_2} -o HLAtyping -n ${PREFIX}
#python ${iTuNES_BIN_PATH}/HLA_extract.py HLAtyping/${PREFIX}_HLAminer_HPTASR.csv HLAtyping/${PREFIX}_HLA_str.txt
#hla_str=`cat HLAtyping/${PREFIX}_HLA_str.txt`
razers3 -i 95 -tc 8 -m 1 -dr 0 -o HLAtyping/finished_filter.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq
samtools bam2fq HLAtyping/finished_filter.bam > HLAtyping/finished_filter.fastq 
rm HLAtyping/finished_filter.bam
python /home/zhouchi/software/OptiType/OptiTypePipeline.py -i HLAtyping/finished_filter.fastq --dna -o HLAtyping/
rm HLAtyping/finished_filter.fastq
if [ -f HLAtyping/${PREFIX}_optitype_hla_type ];then
	rm HLAtyping/${PREFIX}_optitype_hla_type
fi
dir=`ls HLAtyping/`
python ${iTuNES_BIN_PATH}/optitype_ext.py -i HLAtyping/${dir}/${dir}_result.tsv -o HLAtyping -s ${PREFIX}
!
#######VEP Annotation########
:<<!
echo -e "\033[40;32mVEP Annotation......\t"$DATE"\033[0m"

#VEP_DATA=/home/zhouchi/database/Annotation/vep_data
#vep.pl -i somatic_mutation/${PREFIX}/${PREFIX}_SNVs_only.recode.vcf --cache --dir $VEP_DATA --dir_cache #$VEP_DATA --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_vep_ann.txt --force_overwrite
########Preprocess########
#python ${iTuNES_BIN_PATH}/animo_acid_prepare.py -i somatic_mutation/${PREFIX}_vep_ann.txt -o netmhc -s ${PREFIX}
########netMHC###########
#netMHCpan -a $hla_str -f netmhc/${PREFIX}_netmhc_input.fasta -l 8,9,10,11 > netmhc/${PREFIX}_netmhc_result.txt
#######filtering######
set -x
#python ${iTuNES_BIN_PATH}/netMHC_result_parse.py -i netmhc/${PREFIX}_netmhc_result.txt -o netmhc -s ${PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
set +x
#######netchop########
set -x
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_final_neo_candidate.txt -o netctl -s ${PREFIX}
set +x
###indel identification###
#vep.pl -i somatic_mutation/${PREFIX}/${PREFIX}_INDELs_only.recode.vcf --cache --dir $VEP_DATA --dir_cache #$VEP_DATA --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_vep_ann.txt --force_overwrite

!
########Preprocess########
#python ${iTuNES_BIN_PATH}/deletion2fasta.py -i somatic_mutation/${PREFIX}_D_vep_ann.txt -c somatic_mutation/${PREFIX}_D.vcf -o netmhc -s ${PREFIX}
#python ${iTuNES_BIN_PATH}/insertion2fasta.py -i somatic_mutation/${PREFIX}_SI_vep_ann.txt -c somatic_mutation/${PREFIX}_SI.vcf -o netmhc -s ${PREFIX}
########netMHC###########
#netMHCpan -a $hla_str -f netmhc/${PREFIX}_del.fasta -l 9 > netmhc/${PREFIX}_del_netmhc_result.txt
#netMHCpan -a $hla_str -f netmhc/${PREFIX}_ins.fasta -l 9 > netmhc/${PREFIX}_ins_netmhc_result.txt
#######filtering######
#python ${iTuNES_BIN_PATH}/netMHC_result_parse.py -i netmhc/${PREFIX}_del_netmhc_result.txt -o netmhc -s ${DEL_PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
#python ${iTuNES_BIN_PATH}/netMHC_result_parse.py -i netmhc/${PREFIX}_ins_netmhc_result.txt -o netmhc -s ${INS_PREFIX} -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
#######netchop########
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${DEL_PREFIX}_final_neo_candidate.txt -o netctl -s ${DEL_PREFIX}
#python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${INS_PREFIX}_final_neo_candidate.txt -o netctl -s ${INS_PREFIX}


###finished
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"
