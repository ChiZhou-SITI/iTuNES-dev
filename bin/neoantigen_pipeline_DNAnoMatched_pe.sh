#!/bin/bash

#===============================================================================
#
#          FILE:  neoantigen_pipeline_DNAnoMatched_pe.sh
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
	echo "	 This is an pipeline calling mutation from  matched tumor/normal WGSS raw data"
	echo "   Usage:"
	echo "	 bash neoantigen_pipeline_DNAnoMatched_pe.sh [option] <-p tumor.fq(fastq)> <-q tumor.fq(fastq)> <-o out_dir> <-g genome> <-P prefix> <-I bwa_index> <-r reference> <-c cpu_num> <-C coverage> <-b bindingaffinity_cutoff> <-E expression_file> <-f fpkm_cutoff>"
	echo ""
	echo "   Required argument:"
	echo "  	-p tumor fastq1. first fastq file from tumor sample "
	echo "  	-q tumor fastq2. second fastq file from tumor sample "
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

TEMP=`getopt -o p:q:o:g:P:I:r:c:C:a:b:E:f: \
--long T_fq_1:,T_fq_2:,outpath:,genome:,prefix:,bwa_index:,reference:,CPU:Coverage:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff: \
	-n 'neoantigen_pipeline_DNAMatched.sh' -- "$@"`

if [ $? != 0 ];then
	echo "Terminating..." >&2 ; exit 1 ; fi


eval set -- "$TEMP"
while true
do
	case "$1" in
		-p | --T_fq_1) T_fastq_1=$2;shift 2;;
		-q | --T_fq_2) T_fastq_2=$2;shift 2;;
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

#if [ ! -f "T_fastq_1" -o ! -f "$T_fastq_2" ];then
#	echo -e "\033[40;31;1mERROR: please use -p and -q to specify the right WGS fastq file\033[0m"
#	exit 1
#fi
#if [ ! -n "$GENOME" ];then
#	echo -e "\033[40;31;1mERROR: please use -g to specify genome used, eg. mm10\033[0m"
#	exit 1
#fi

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
ref_fasta='/home/zhouchi/database/hg19_ref/hg19.fasta'
SortSam="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard SortSam"
ref_prefix='/home/zhouchi/database/hg19_ref/hg19'
echo -e "\033[40;32mcalling somatic mutaiton using GATK......\t"$DATE"\033[0m"
set -x

BWA_pefq2bam(){
    local input_fq1=$1
    local input_fq2=$2
    local output_bam=$3
    local output_dup=$4
    sam=`${iTuNES_BIN_PATH}/newtf`.sam;	#mkfifo $sam;
    bwa mem -t 10 $ref_prefix $input_fq1 $input_fq2 > $sam                      # Make the output of bwa as bam
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
    $RealignerTargetCreator -R $ref_fasta -I $input_bam -known $known_indels1 -known $known_indels2 -filterRNC -o $intervals 
    $IndelRealigner         -R $ref_fasta -I $input_bam -known $known_indels1 -known $known_indels2 -targetIntervals $intervals -filterRNC -o $output      # The parallel -nt -nct are not supported
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
    $BaseRecalibrator  -R $ref_fasta -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3 -o $before_table
    $PrintReads        -R $ref_fasta -I $input_bam -BQSR $before_table -o $output
    #$BaseRecalibrator  -R $ref_fasta -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3  -BQSR $before_table -o $after_table
    #$AnalyzeCovariates -R $ref_fasta -before $before_table -after $after_table -plots $recal_plots -l DEBUG        
    rm $before_table
    # note that Rscript will be executed. many packages are needed. add -l DEBUG to check which library are needed 
    # if the execution failed. The following library are required:
    # reshape, gplots, ggplot2, bitops, caTools, colorspace, gtools, gData/Storage1, RColorBrewer, gsalib
}



callvariants_dna(){
    local input_dna_bam=$1
    local output=$2
    $HaplotypeCaller -R $ref_fasta -I $input_dna_bam -o $output -stand_call_conf 30 --dbsnp $dbsnp138
}


##VQSR
variantrecalibration(){
    local input_vcf=$1
    local output_vcf=$2
    local output_snp_recal_R=$3
    local output_indel_recal_R=$4
    local snpTranchesPdf=$5
    local indelTranchesPdf=$6

    local recal=`${iTuNES_BIN_PATH}/newtf`.recal
    local tranches=`${iTuNES_BIN_PATH}/newtf`.tranches
    $VariantRecalibrator -R $ref_fasta -input $input_vcf \
        -recalFile    $recal \
        -tranchesFile $tranches \
        -rscriptFile  $output_snp_recal_R \
        -mode SNP \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
        -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $OneKG \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp138 \
        -an MQRankSum -an ReadPosRankSum -an FS -an DP
        #--minNumBadVariants 100

    local recalSNPVCF=`${iTuNES_BIN_PATH}/newtf`.vcf
    $ApplyRecalibration  -R $ref_fasta -input $input_vcf \
        -recalFile    $recal \
        -tranchesFile $tranches \
        -o            $recalSNPVCF \
        -mode SNP \
        -ts_filter_level 99.0 

    mv $tranches.pdf $snpTranchesPdf

    $VariantRecalibrator -R $ref_fasta -input $recalSNPVCF\
        -recalFile    $recal \
        -tranchesFile $tranches \
        -rscriptFile  $output_indel_recal_R \
        -mode INDEL \
        --maxGaussians 2  \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp138 \
        -an DP -an FS -an ReadPosRankSum -an MQRankSum 
        #-minNumBadVariants 100

    $ApplyRecalibration  -R $ref_fasta -input $recalSNPVCF \
        -recalFile    $recal \
        -tranchesFile $tranches \
        -o            $output_vcf \
        -mode INDEL \
        -ts_filter_level 99.0

    mv $tranches.pdf $indelTranchesPdf

    rm $recal* $tranches* $recalSNPVCF*
}
dbsnp138='/home/zhouchi/database/Annotation/hg38_vcf/dbsnp_138.hg38.vcf'
hapmap='/home/zhouchi/database/Annotation/hg38_vcf/hapmap_3.3.hg38.vcf'
omni='/home/zhouchi/database/Annotation/hg38_vcf/1000G_phase1.snps.high_confidence.hg38.vcf'
OneKG='/home/zhouchi/database/Annotation/hg38_vcf/1000G_omni2.5.hg38.vcf'
mills='/home/zhouchi/database/Annotation/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf'
#BWA_pefq2bam	$T_fastq_1	$T_fastq_2	alignments/${PREFIX}_tumor_bwa.bam	alignments/${PREFIX}_tumor_bwa_dedup.dup
#realignment             alignments/${PREFIX}_tumor_bwa.bam        $OneKG     $mills  alignments/${PREFIX}_tumor_bwa_realign.bam
#baserecalibration       alignments/${PREFIX}_tumor_bwa_realign.bam      $dbsnp138      $OneKG     $mills   alignments/${PREFIX}_tumor_bwa_realign_recal.bam
#callvariants_dna            alignments/${PREFIX}_tumor_bwa_realign_recal.bam	somatic_mutation/${PREFIX}_raw.vcf
#variantrecalibration somatic_mutation/${PREFIX}_raw.vcf somatic_mutation/${PREFIX}_VQSR.vcf somatic_mutation/${PREFIX}_snp_recal.R somatic_mutation/${PREFIX}_indel_recal.R somatic_mutation/${PREFIX}_snpTranches.pdf somatic_mutation/${PREFIX}_indelTranches.pdf
set +x


##########neoantigen identification#####
#######HLA TYPING#######
echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
#python $HLAtyping -i ${T_fastq_1} ${T_fastq_2} -o HLAtyping -n ${PREFIX}
#python ${iTuNES_BIN_PATH}/HLA_extract.py HLAtyping/${PREFIX}_HLAminer_HPTASR.csv HLAtyping/${PREFIX}_HLA_str.txt
#hla_str=`cat HLAtyping/${PREFIX}_HLA_str.txt`
#cd HLAtyping
#python /home/zhouchi/software/seq2HLA/seq2HLA.py -1 ${T_fastq_1} -2 ${T_fastq_2} -r ${PREFIX}seq2HLA -p 4
#python ${iTuNES_BIN_PATH}/seq2HLAExt.py -i ${PREFIX}seq2HLA-ClassI.HLAgenotype4digits -o . -s ${PREFIX}
#rm ${PREFIX}seq2HLA*
#cd ..
#hla_str=`cat HLAtyping/${PREFIX}_hlatype_result`
#######VEP Annotation########
echo -e "\033[40;32mVEP Annotation......\t"$DATE"\033[0m"

#VEP_DATA=/data/PUBLIC/vep_data
#vep.pl -i somatic_mutation/${PREFIX}_VQSR.vcf --cache --dir $VEP_DATA --dir_cache #$VEP_DATA --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_vep_ann.txt --force_overwrite
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
echo -e "\033[40;32mindel calling ......\t"$DATE"\033[0m"
samtools index alignments/${PREFIX}_normal_mkdup_filter.bam
samtools index alignments/${PREFIX}_tumor_mkdup_filter.bam
cd somatic_mutation
DEL_PREFIX=${PREFIX}_del
INS_PREFIX=${PREFIX}_ins
echo -e ../alignments/${PREFIX}_tumor_mkdup_filter.bam'\t250\t'${PREFIX} > ${PREFIX}_config.txt
echo -e ../alignments/${PREFIX}_normal_mkdup_filter.bam'\t250\t'${PREFIX} >> ${PREFIX}_config.txt
pindel -f /home/zhouchi/database/Annotation/Fasta/human.fasta -i ${PREFIX}_config.txt -c ALL -o ${PREFIX}
pindel2vcf -p ${PREFIX}_D -r /home/zhouchi/database/Annotation/Fasta/human.fasta -R ucsc_hg19 -d 20090201 -v ${PREFIX}_D.vcf
vep.pl -i ${PREFIX}_D.vcf --cache --dir $VEP_DATA --dir_cache $VEP_DATA --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is coding_sequence_variant" -o ${PREFIX}_D_vep_ann.txt --force_overwrite
pindel2vcf -p ${PREFIX}_SI -r /home/zhouchi/database/Annotation/Fasta/human.fasta -R ucsc_hg19 -d 20090201 -v ${PREFIX}_SI.vcf
vep.pl -i ${PREFIX}_SI.vcf --cache --dir $VEP_DATA --dir_cache $VEP_DATA --force_overwrite  --symbol -o STDOUT --offline | filter_vep.pl --ontology --filter "Consequence is coding_sequence_variant" -o ${PREFIX}_SI_vep_ann.txt --force_overwrite
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mfinish indel calling and annotation ......\t"$DATE"\033[0m"
cd ..
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
