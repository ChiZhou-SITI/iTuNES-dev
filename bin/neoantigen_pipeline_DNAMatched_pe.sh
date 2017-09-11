#!/bin/bash

#===============================================================================
#
#          FILE:  neoantigen_pipeline_DNAMatched.sh
#
# 
#   DESCRIPTION:  This is an pipeline calling mutation from  matched tumor/normal WGSS raw data
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
set -e
#######--Argument--######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: ChiZhou"
	echo "	 This is an pipeline calling mutation from  matched tumor/normal WGSS raw data"
	echo "   Usage:"
	echo "	 bash neoantigen_pipeline_DNAMatched.sh [option] <-p normal.fq(fastq)> <-s tumor.fq(fastq)> <-g genome> <-q normal_2.fq> <-t tumor_2.fq> <-o out_dir> <-g genome> <-P prefix> <-I bwa_index> <-r reference> <-c cpu_num> <-C coverage> <-b bindingaffinity_cutoff> <-E expression_file> <-f >"
	echo ""
	echo "   Required argument:"
	echo "  	-p normal fastq1.first fastq file from normal sample "
	echo "  	-q normal fastq2.second fastq file from normal sample "
	echo "  	-s tumor fastq1.first fastq file from tumor sample "
	echo "  	-t tumor fastq2.first fastq file from tumor sample "
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
	echo -e "	The program can deal with pair-end matched tumor/normal wgs data. If not, please modify the program."
	echo "	Program needed: bwa edwstates picards samtools bedtools Varscan VEP HLAminer netMHC netchop Sequenza PyClone "
	echo "	Make sure all the program is contained by environment viariable PATH."
	echo "	The pipeline also need index files builded already, if not the pipeline will build it."
	echo ""
	echo "    (-^V^-) Enjoy yourself~~~"
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo ""
}
if [ $# -lt 2 ];then
	help_info
	exit 1
fi
TEMP=`getopt -o p:q:s:t:o:g:P:I:r:c:C:a:b:E:f:l: \
--long N_fq_1:,N_fq_2:,T_fq_1:,T_fq_2:,outpath:,genome:,prefix:,bwa_index:,reference:,CPU:Coverage:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff:,hla_type: \
	-n 'neoantigen_pipeline_DNAMatched.sh' -- "$@"`
if [ $? != 0 ];then
	echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
while true
do
	case "$1" in
		-p | --N_fq_1) N_fastq_1=$2;shift 2;;
		-q | --N_fq_2) N_fastq_2=$2;shift 2;;
		-s | --T_fq_1) T_fastq_1=$2;shift 2;;
		-t | --T_fq_2) T_fastq_2=$2;shift 2;;
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
		-l | --hla_type) Hla_Type=$2;shift 2;;
		--) shift;break;;
		*) echo "Internal error!";exit 1;;
	esac

done
###########check Argument setting###########
#if [ ! -f "$N_fastq1" -o ! -f "T_fastq1" -o ! -f "$N_fastq2" -o ! -f "$T_fastq2" ];then
#	echo -e "\033[40;31;1mERROR: please use -p.-q,-s and -t to specify the right WGS fastq file\033[0m"
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
if [ "$BWA_INDEX" = "None" ];then
	BWA_INDEX=/home/zhouchi/database/Annotation/Index/bwa/${GENOME}
fi
if [ "$REFERENCE" = "None" ];then
	REFERENCE=/home/zhouchi/database/Annotation/Fasta/${GENOME}.fasta
fi
if [ ! -n "$COVERAGE" ];then
	COVERAGE=0
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
###check index
if [ -f ${BWA_INDEX}/${GENOME}".amb" -a -f ${BWA_INDEX}/${GENOME}".ann" -a -f ${BWA_INDEX}/${GENOME}".bwt" -a -f ${BWA_INDEX}/${GENOME}".pac" -a -f ${BWA_INDEX}/${GENOME}".sa" ];then
	echo -e "\033[40;35;1mIndex: "${BWA_INDEX}"\033[0m"
else
	echo "no index file, make bwa index"
	if [ -f ${GENOME_FA} ];then
		DATE=`date --date="-24 hour"`
		echo -e "\033[40;32mmake index......\t"$DATE"\033[0m"
		set -x
		mkdir -p ${BWA_INDEX}
		bwa index -p ${BWA_INDEX}/${GENOME} -a bwtsw ${REFERENCE}
		set +x
	else
		echo -e "\033[40;31;1mERROR: no index file and no genome.fa to build it\033[0m"
		exit 1
	fi
	echo -e "\033[40;35;1mIndex: "${BWA_INDEX}"\033[0m"
fi

Picard=/home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar
GATK="/home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar"
BuildBamIndex="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard BuildBamIndex"
MarkDuplicates="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard MarkDuplicates"
AddOrReplaceReadGroups="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $Picard AddOrReplaceReadGroups"
RealignerTargetCreator="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T RealignerTargetCreator -nt 8 -dt NONE"
IndelRealigner="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T IndelRealigner"
BaseRecalibrator="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T BaseRecalibrator -nct 8"
AnalyzeCovariates="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T AnalyzeCovariates -dt NONE"
HaplotypeCaller="java -Djava.io.tmpdir=/home/zhouchi/tmp -Xmx8g -jar $GATK -T HaplotypeCaller -nct 8 -dt NONE"
SplitNCigarReads="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T SplitNCigarReads"
PrintReads="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T PrintReads -nct 8 -dt NONE"
VariantFiltration="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T VariantFiltration"
GATK_Mutect2="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T MuTect2 -nct 8"
SelectVariants="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T SelectVariants"
####################
###Align with bwa###
####################
###align

dbsnp138='/home/zhouchi/database/Annotation/hg38_vcf/dbsnp_138.hg38.vcf'
hapmap='/home/zhouchi/database/Annotation/hg38_vcf/hapmap_3.3.hg38.vcf'
omni='/home/zhouchi/database/Annotation/hg38_vcf/1000G_phase1.snps.high_confidence.hg38.vcf'
OneKG='/home/zhouchi/database/Annotation/hg38_vcf/1000G_omni2.5.hg38.vcf'
mills='/home/zhouchi/database/Annotation/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf'
cosmic='/home/zhouchi/database/Annotation/hg38_vcf/CosmicCodingMuts_chr_M_sorted.vcf'
######check other files########
hlatyping(){
	local PREFIX=$1
	rm -rf HLAtyping/*
	#razers3 -i 95 -tc $CPU -m 1 -dr 0 -o HLAtyping/fished_test_1.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_1
	#razers3 -i 95 -tc $CPU -m 1 -dr 0 -o HLAtyping/fished_test_2.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_2
	#samtools bam2fq HLAtyping/fished_test_1.bam > HLAtyping/test_1_fished.fastq 
	#samtools bam2fq HLAtyping/fished_test_2.bam > HLAtyping/test_2_fished.fastq
	#rm HLAtyping/fished_test_1.bam HLAtyping/fished_test_2.bam
	python /home/zhouchi/software/OptiType/OptiTypePipeline.py -i $T_fastq_1 $T_fastq_2 --dna -o HLAtyping/
	#rm HLAtyping/test_1_fished.fastq HLAtyping/test_2_fished.fastq
	dir=`ls HLAtyping/`
	python ${iTuNES_BIN_PATH}/optitype_ext.py -i HLAtyping/${dir}/${dir}_result.tsv -o HLAtyping -s ${PREFIX}
}
mapping_qc_normal(){
	local $PREFIX
	bwa mem -t ${CPU} ${BWA_INDEX}/${GENOME} ${N_fastq_1} ${N_fastq_2} > alignments/tmp_${PREFIX}_normal.sam

	samtools view -bhS -@ ${CPU} alignments/tmp_${PREFIX}_normal.sam > alignments/tmp_${PREFIX}_normal.bam 

	samtools sort -@ ${CPU} -m 8G alignments/tmp_${PREFIX}_normal.bam alignments/${PREFIX}_normal_unfilter 

	
	samtools index alignments/${PREFIX}_normal_unfilter.bam 
	samtools view -b -@ ${CPU} alignments/${PREFIX}_normal_unfilter.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > alignments/${PREFIX}_normal.bam 
	samtools index alignments/${PREFIX}_normal.bam 
	set +x
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mrun picard mark duplicates on non-UMI......\t"$DATE"\033[0m"
	set -x
	time java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar MarkDuplicates \
		INPUT=alignments/${PREFIX}_normal.bam OUTPUT=alignments/${PREFIX}_normal_mkdup.bam \
		METRICS_FILE=${PREFIX}_normal_dup_qc.txt ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	set +x

	###fliter bam on flags and threashold
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mfilter on flag 512 and threadshold score: "$SCORE"......\t"$DATE"\033[0m"
	samtools view -F 512 -b -@ ${CPU} alignments/${PREFIX}_normal_mkdup.bam > alignments/${PREFIX}_normal_mkdup_filter.bam
	samtools index alignments/${PREFIX}_normal_mkdup_filter.bam
	rm alignments/${PREFIX}_normal_mkdup.bam
	rm alignments/tmp_${PREFIX}_normal.bam alignments/tmp_${PREFIX}_normal.sam
	rm alignments/${PREFIX}_normal_unfilter.bam alignments/${PREFIX}_normal_unfilter.bam.bai 

}
mapping_qc_tumor(){
	local $PREFIX

	bwa mem -t ${CPU} ${BWA_INDEX}/${GENOME} ${T_fastq_1} ${T_fastq_2} > alignments/tmp_${PREFIX}_tumor.sam 

	samtools view -bhS -@ ${CPU} alignments/tmp_${PREFIX}_tumor.sam > alignments/tmp_${PREFIX}_tumor.bam

	samtools sort -@ ${CPU} -m 8G alignments/tmp_${PREFIX}_tumor.bam alignments/${PREFIX}_tumor_unfilter

	
	samtools index alignments/${PREFIX}_tumor_unfilter.bam 

	samtools view -b -@ ${CPU} alignments/${PREFIX}_tumor_unfilter.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > alignments/${PREFIX}_tumor.bam
	samtools index alignments/${PREFIX}_tumor.bam
	
	set +x

	#######################
	###Filter alignments###
	#######################
	###mark duplicates with picard
	DATE=`date --date="-24 hour"`

	echo -e "\033[40;32mrun picard mark duplicates on non-UMI......\t"$DATE"\033[0m"
	set -x
	time java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar MarkDuplicates \
		INPUT=alignments/${PREFIX}_tumor.bam OUTPUT=alignments/${PREFIX}_tumor_mkdup.bam \
		METRICS_FILE=${PREFIX}_tumor_dup_qc.txt ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	
	set +x

	###fliter bam on flags and threashold
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mfilter on flag 512 and threadshold score: "$SCORE"......\t"$DATE"\033[0m"
	samtools view -F 512 -b -@ ${CPU} alignments/${PREFIX}_tumor_mkdup.bam > alignments/${PREFIX}_tumor_mkdup_filter.bam
	samtools index alignments/${PREFIX}_tumor_mkdup_filter.bam
	rm alignments/${PREFIX}_tumor_mkdup.bam
	rm alignments/tmp_${PREFIX}_tumor.bam alignments/tmp_${PREFIX}_tumor.sam
	rm  alignments/${PREFIX}_tumor_unfilter.bam alignments/${PREFIX}_tumor_unfilter.bam.bai
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
    local before_table=`${iTuNES_BIN_PATH}/newtf`
    $BaseRecalibrator  -R $REFERENCE -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3 -o alignments/$before_table
    $PrintReads        -R $REFERENCE -I $input_bam -BQSR alignments/$before_table -o $output       
    rm alignments/$before_table
}
GATK_mutation_calling(){
	local tumor_bam=$1
	local normal_bam=$2
	local dbsnp138=$3
	local cosmic=$4
	local output_vcf=$5
	$GATK_Mutect2 \
     -R $REFERENCE \
     -I:tumor ${tumor_bam} \
     -I:normal ${normal_bam} \
     --dbsnp $dbsnp138 \
     --cosmic $cosmic \
     -o ${output_vcf}
}
GATK_pipe(){
	local PREFIX=$1
	$AddOrReplaceReadGroups I=alignments/${PREFIX}_normal_mkdup_filter.bam O=alignments/${PREFIX}_normal_mkdup_filter_add.bam SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20      # Add your own annotation
	$AddOrReplaceReadGroups I=alignments/${PREFIX}_tumor_mkdup_filter.bam O=alignments/${PREFIX}_tumor_mkdup_filter_add.bam SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20      # Add your own annotation
	$BuildBamIndex I=alignments/${PREFIX}_normal_mkdup_filter_add.bam O=alignments/${PREFIX}_normal_mkdup_filter_add.bam.bai VALIDATION_STRINGENCY=SILENT
	$BuildBamIndex I=alignments/${PREFIX}_tumor_mkdup_filter_add.bam O=alignments/${PREFIX}_tumor_mkdup_filter_add.bam.bai VALIDATION_STRINGENCY=SILENT
	realignment alignments/${PREFIX}_normal_mkdup_filter_add.bam $OneKG $mills alignments/${PREFIX}_normal_mkdup_filter_add_realign.bam
	realignment alignments/${PREFIX}_tumor_mkdup_filter_add.bam $OneKG $mills alignments/${PREFIX}_tumor_mkdup_filter_add_realign.bam
	baserecalibration alignments/${PREFIX}_normal_mkdup_filter_add_realign.bam $dbsnp138 $OneKG $mills alignments/${PREFIX}_normal_recal.bam
	baserecalibration alignments/${PREFIX}_tumor_mkdup_filter_add_realign.bam $dbsnp138 $OneKG $mills alignments/${PREFIX}_tumor_recal.bam
	#GATK_mutation_calling alignments/${PREFIX}_tumor_recal.bam alignments/${PREFIX}_normal_recal.bam $dbsnp138 $cosmic somatic_mutation/${PREFIX}_gatk_vep_input.vcf
}
varscan_somatic_calling(){
	local PREFIX=$1
	rm -rf somatic_mutation/*
	cd somatic_mutation
	mkfifo ${PREFIX}_normal.fifo
	mkfifo ${PREFIX}_tumor.fifo
	samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ../alignments/${PREFIX}_normal_recal.bam  > ${PREFIX}_normal.fifo &
	samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ../alignments/${PREFIX}_tumor_recal.bam  > ${PREFIX}_tumor.fifo &
	java -jar /home/zhouchi/software/varscan/VarScan.v2.4.2.jar somatic ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX} #--output-vcf 1
	java -jar /home/zhouchi/software/varscan/VarScan.v2.4.2.jar processSomatic ${PREFIX}.snp
	rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
	cd ..
}


varscansnv2fasta(){
	local VEP_CACHE=$1
	local PREFIX=$2
	sed '1d' somatic_mutation/${PREFIX}.snp.Somatic | awk -F '\t' '{print $1,$2,$2,$3"/"$4}' > somatic_mutation/${PREFIX}_snv_vep_input.vcf
	vep -i somatic_mutation/${PREFIX}_snv_vep_input.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol --offline -o somatic_mutation/${PREFIX}_snv_all_vep_ann.txt
	vep -i somatic_mutation/${PREFIX}_snv_vep_input.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_snv_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/snv2fasta.py -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -o netmhc -s ${PREFIX}

}
varscanindel2fasta(){
	local VEP_CACHE=$1
	local PREFIX=$2
	python ${iTuNES_BIN_PATH}/varscan_indel_preprocess.py -i somatic_mutation/${PREFIX}.indel -o somatic_mutation -s ${PREFIX}
	vep -i somatic_mutation/${PREFIX}_varscan_indel.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt -o netmhc -s ${PREFIX}_varscan
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_varscan
}
strelka_indel_calling(){
	local $PREFIX
	if [ -d strelka_indel ];then
		rm -rf strelka_indel
	fi
	python /usr/local/bin/configureStrelkaSomaticWorkflow.py --tumorBam=alignments/${PREFIX}_tumor_recal.bam --normalBam=alignments/${PREFIX}_normal_recal.bam --referenceFasta=${REFERENCE} --config=/usr/local/bin/configureStrelkaSomaticWorkflow.py.ini --runDir=strelka_indel --exome
	python strelka_indel/runWorkflow.py -m local -j $CPU -q ${PREFIX}_strelka -g 32 --quiet
}
strelkaindel2fasta(){
	local VEP_CACHE=$1
	local PREFIX=$2
	vep -i strelka_indel/results/variants/somatic.indels.vcf.gz --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt -o netmhc -s ${PREFIX}_strelka
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_strelka
}

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

snv_fasta2netmhcpan2netCTLPAN(){
	local PREFIX=$1
	local Exp_file=$2
	local Binding_Aff_Fc_Cutoff=$3
	local Binding_Aff_Cutoff=$4 
	local Fpkm_Cutoff=$5
	local Hla_Type=$6
	netmhcpan netmhc/${PREFIX}_snv.fasta ${hla_str} ${PREFIX}_snv_netmhc.txt netmhc 100
	python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_snv_netmhc.txt -g netmhc/${PREFIX}_snv.fasta -o netmhc -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
	python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_snv_final_neo_candidate.txt -o netctl -s ${PREFIX}_snv

}
strelka_varscanindel_fasta2netmhcpan2netCTLPAN(){
	local PREFIX=$1
	local Exp_file=$2
	local Binding_Aff_Fc_Cutoff=$3
	local Binding_Aff_Cutoff=$4 
	local Fpkm_Cutoff=$5
	local Hla_Type=$6
	cat netmhc/${PREFIX}_strelka_del.fasta > netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_strelka_ins.fasta >> netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_varscan_del.fasta >> netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_varscan_ins.fasta >> netmhc/${PREFIX}_indel.fasta
	netmhcpan netmhc/${PREFIX}_indel.fasta ${hla_str} ${PREFIX}_indel_netmhc.txt netmhc 100
	python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_indel_netmhc.txt -g netmhc/${PREFIX}_indel.fasta -o netmhc -s ${PREFIX}_indel -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
	python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_indel_final_neo_candidate.txt -o netctl -s ${PREFIX}_indel
	}

pyclone_copynumber_varscan(){
	local PREFIX=$1
	chrlist=`sed '1d' somatic_mutation/${PREFIX}.snp | cut -f1 | uniq`
	echo 'chromosome	position	base.ref	depth.normal	depth.tumor	depth.ratio	Af	Bf	zygosity.normal	GC.percent	good.reads	AB.normal	AB.tumor	tumor.strand' > pyclone/${PREFIX}.seqz
	for s in $chrlist
	do
	{
		str_chr='-C '$s
		python ${iTuNES_BIN_PATH}/sequenza-utils.py bam2seqz -gc /home/zhouchi/database/Annotation/human.gc50Base.txt.gz $str_chr --fasta $REFERENCE -n alignments/${PREFIX}_normal_recal.bam -t alignments/${PREFIX}_tumor_recal.bam | sed '1d' >> pyclone/${s}.seqz
	} &
	done
	wait
	
	for s in $chrlist
	do
	{
		cat pyclone/${s}.seqz >> pyclone/${PREFIX}.seqz
		rm pyclone/${s}.seqz
	}
	done
	gzip -f pyclone/${PREFIX}.seqz > pyclone/${PREFIX}.seqz.gz
	python ${iTuNES_BIN_PATH}/sequenza-utils.py seqz-binning -w 50 -s pyclone/${PREFIX}.seqz.gz | gzip > pyclone/${PREFIX}.small.seqz.gz
	Rscript ${iTuNES_BIN_PATH}/sequenza_test.R pyclone/${PREFIX}.small.seqz.gz pyclone/ ${PREFIX}
	python ${iTuNES_BIN_PATH}/pyclone_input.py -n netctl/${PREFIX}_snv_netctl_concact.txt -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -s somatic_mutation/${PREFIX}.snp.Somatic -c pyclone/${PREFIX}_seg_copynumber.txt -o pyclone -S ${PREFIX} -C ${COVERAGE}
	TUMOR_CONTENT=`cat pyclone/${PREFIX}_cellularity.txt`
	PyClone run_analysis_pipeline --in_files pyclone/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir pyclone 
	python ${iTuNES_BIN_PATH}/neo_pyclone_annotation.py -n netctl/${PREFIX}_snv_netctl_concact.txt -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -s pyclone/tables/loci.tsv -o netctl -S ${PREFIX}
}
pipeline_all(){
	set -x
	local PREFIX=$1
	if [ "$Hla_Type" = 'None' ];then
		DATE=`date --date="-24 hour"`
		echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
		begin_time=$(date +%s)
		#hlatyping $PREFIX
		end_time=$(date +%s)
		cost_time=$(($begin_time - $end_time))
		echo "it spended $cost_time seconds running the hlatyping." >> log_file/running_time.txt
		#hla_str=`cat HLAtyping/${PREFIX}_optitype_hla_type`
	else
		hla_str=${Hla_Type}
	fi
	
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mmapping......\t"$DATE"\033[0m"
	begin_time=$(date +%s)
	mapping_qc_normal $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds running the normal sample mapping and quality control." >> log_file/running_time.txt
	begin_time=$(date +%s)
	mapping_qc_tumor $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds running the tumor sample mapping and quality control." >> log_file/running_time.txt

	begin_time=$(date +%s)
	GATK_pipe $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds running the GATK process." >> log_file/running_time.txt
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mcalling somatic mutaiton and copynumber profile......\t"$DATE"\033[0m"
	begin_time=$(date +%s)
	varscan_somatic_calling $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds to call somatic profile." >> log_file/running_time.txt
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32mindel calling using stelka......\t"$DATE"\033[0m"

	begin_time=$(date +%s)
	strelka_indel_calling $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds to use stelka to call indel." >> log_file/running_time.txt
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32fasta generate ......\t"$DATE"\033[0m"
	VEP_CACHE=/home/zhouchi/database/Annotation/vep_data/
	varscansnv2fasta /home/zhouchi/database/Annotation/vep_data $PREFIX
	varscanindel2fasta /home/zhouchi/database/Annotation/vep_data $PREFIX
	strelkaindel2fasta /home/zhouchi/database/Annotation/vep_data $PREFIX
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32Runnning netmhcpan and netCTLPAN ......\t"$DATE"\033[0m"
	####snv
	begin_time=$(date +%s)
	snv_fasta2netmhcpan2netCTLPAN $PREFIX ${Exp_file} ${Binding_Aff_Fc_Cutoff} ${Binding_Aff_Cutoff} ${Fpkm_Cutoff} ${Hla_Type}
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds to run netMHCpan for snv." >> log_file/running_time.txt
	####indel
	begin_time=$(date +%s)
	strelka_varscanindel_fasta2netmhcpan2netCTLPAN $PREFIX ${Exp_file} ${Binding_Aff_Fc_Cutoff} ${Binding_Aff_Cutoff} ${Fpkm_Cutoff} ${Hla_Type}
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds to run netMHCpan for indel." >> log_file/running_time.txt
	DATE=`date --date="-24 hour"`
	echo -e "\033[40;32Runnning pyclone to get cellularity preverence......\t"$DATE"\033[0m"
	begin_time=$(date +%s)
	pyclone_copynumber_varscan $PREFIX
	end_time=$(date +%s)
	cost_time=$(($begin_time - $end_time))
	echo "it spended $cost_time seconds to run pyclone." >> log_file/running_time.txt

	set +x
}
pipeline_all $PREFIX
DATE=`date --date="-24 hour"`
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"