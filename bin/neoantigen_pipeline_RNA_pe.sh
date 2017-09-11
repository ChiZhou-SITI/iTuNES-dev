
#!/bin/bash

#===============================================================================
#
#          FILE:  neoantigen_pipeline_RNAnoMatched_pe.sh
#
# 
#   DESCRIPTION:  This is an pipeline calling mutation from tumor RNA-seq raw data
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   ChiZhou 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  3/31/2017 15:09 CST
#      REVISION:  ---
#===============================================================================
set -e
#######--Argument--######
help_info(){
	echo ""
	echo -e "\033[32m  =========================================================================================================================\033[0m"
	echo "   Author: ChiZhou"
	echo "	 This is an pipeline calling mutation from  tumor RNA-seq raw data"
	echo "   Usage:"
	echo "	 bash neoantigen_pipeline_RNA_pe.sh [option] <-p tumor.fq(fastq)> <-q tumor.fq(fastq)> <-o out_dir> <-g genome> <-P prefix> <gtf_file> <-I star_index> <-r reference> <-c cpu_num> <-C coverage> <-b bindingaffinity_cutoff> <-E expression_file> <-f fpkm_cutoff> <-F genefusion_cutoff>"
	echo ""
	echo "   Required argument:"
	echo "  	-p tumor fastq1. first fastq file from tumor sample "
	echo "  	-q tumor fastq2. second fastq file from tumor sample "
	echo "  	-g genome species.eg :human "
	
	echo "   Optional arguments:"
	echo "  	-I  STAR index. default: /home/zhouchi/database/Annotation/Index/STAR/genome"
	echo "  	-P prefix for output file. default: sample name"
	echo "  	-o output directory. default: ./"
	echo "  	-r genome reference file. default: /home/zhouchi/database/Annotation/Fasta/genome.fa"
	echo "  	-c CPU thread. --default: 4"
	echo "  	-C site reads coverge cutoff. default: 200"
	echo " 		-a binding affinity fold change cutoff .default: 1"
	echo "  	-b binding affinity cutoff .default: 500"
	echo "  	-E expression profile of tumor sample"
	echo "  	-f gene expression fpkm value cutoff, this option should not set when -E is empty"
	echo "		-F genefusion cutoff for ericscript gene fusion identification"
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

TEMP=`getopt -o p:q:o:g:G:P:I:r:c:C:a:b:E:f:F:l: \
--long T_fq_1:,T_fq_2:,outpath:,genome:,gtf:,prefix:,star_index:,reference:,CPU:Coverage:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff:,genefusion_cutoff:,hla_type: \
	-n 'neoantigen_pipeline_RNA_pe.sh' -- "$@"`

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
		-G | --gtf) GTF=$2;shift 2;;
		-P | --prefix) PREFIX=$2;shift 2;;
		-I | --star_index) STAR_INDEX=$2;shift 2;;
		-r | --reference) REFERENCE=$2;shift 2;;
		-c | --CPU) CPU=$2;shift 2;;
		-C | --Coverage) COVERAGE=$2;shift 2;;
		-b | --binding_aff_cutoff) Binding_Aff_Cutoff=$2;shift 2;;
		-a | --binding_aff_fc_cutoff) Binding_Aff_Fc_Cutoff=$2;shift 2;;
		-E | --expression_file) Exp_file=$2;shift 2;;
		-f | --fpkm_cutoff) Fpkm_Cutoff=$2;shift 2;;
		-F | --genefusion_cutoff) GeneFusion_Cutoff=$2;shift 2;;
		-l | --hla_type) Hla_Type=$2;shift 2;;
		--) shift;break;;
		*) echo "Internal error!";exit 1;;
	esac

done
###########check Argument setting###########
if [ ! -f "$T_fastq_1" -o ! -f "$T_fastq_2" ];then
	echo -e "\033[40;31;1mERROR: please use -p and -q to specify the right WGS fastq file\033[0m"
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

if [ "$STAR_INDEX" = "None" ];then
	STAR_INDEX=/home/zhouchi/database/Annotation/Index/STAR/${GENOME}
fi

if [ "$REFERENCE" = "None" ];then
	REFERENCE=/home/zhouchi/database/Annotation/Fasta/${GENOME}.fasta
fi

if [ "$GTF" = "None" ];then
	GTF=/home/zhouchi/database/Annotation/GTF/${GENOME}.gtf
fi

if [ ! -n "$COVERAGE" ];then
	COVERAGE=100
#elif [ $COVERAGE -lt 0 ];then
#	echo "the COVERAGE should be greater than 0!"
#	exit 1
fi
if [ ! -n "$GeneFusion_Cutoff" ];then
	GeneFusion_Cutoff=0.5
fi
echo $GeneFusion_Cutoff
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

######result fold prepparation#######
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

if [ ! -d expression ];then
	mkdir expression
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

if [ ! -d gene_fusion ];then
	mkdir gene_fusion
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

if [ ! -d annotation ];then
	mkdir annotation
fi



###check index
if [ -f ${STAR_INDEX}/"transcriptInfo.tab" -a -f ${STAR_INDEX}/"sjdbList.out.tab" -a -f ${STAR_INDEX}/"sjdbList.fromGTF.out.tab" -a -f ${STAR_INDEX}/"sjdbInfo.txt" -a -f ${STAR_INDEX}/"SAindex" -a -f ${STAR_INDEX}/"SA" -a -f ${STAR_INDEX}/"genomeParameters.txt" -a -f ${STAR_INDEX}/"Genome" -a -f ${STAR_INDEX}/"geneInfo.tab" -a -f ${STAR_INDEX}/"exonInfo.tab" -a -f ${STAR_INDEX}/"exonGeTrInfo.tab" -a -f ${STAR_INDEX}/"chrStart.txt" -a -f ${STAR_INDEX}/"chrNameLength.txt" -a -f ${STAR_INDEX}/"chrName.txt" -a -f ${STAR_INDEX}/"chrLength.txt" ];then
	echo -e "\033[40;35;1mIndex: "${STAR_INDEX}"\033[0m"
else
	echo "no index file, make STAR index"
	if [ -f ${REFERENCE} ];then
		DATE=`date --date="-24 hour"`
		echo -e "\033[40;32mmake index......\t"$DATE"\033[0m"
		set -x
		STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${STAR_INDEX} --genomeFastaFiles ${REFERENCE} --sjdbGTFfile ${GTF} --sjdbOverhang 99
		set +x
	else
		echo -e "\033[40;31;1mERROR: no index file and no genome.fa to build it\033[0m"
		exit 1
	fi
	echo -e "\033[40;35;1mIndex: "${STAR_INDEX}"\033[0m"
fi


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
SplitNCigarReads="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T SplitNCigarReads"
PrintReads="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T PrintReads -nct 10 -dt NONE"
VariantFiltration="java -Djava.io.tmpdir=/home/zhouchi/tmp -jar $GATK -T VariantFiltration"
####################
###Align with bwa###
####################
###align

dbsnp138='/home/zhouchi/database/Annotation/hg38_vcf/dbsnp_138.hg38.vcf'
hapmap='/home/zhouchi/database/Annotation/hg38_vcf/hapmap_3.3.hg38.vcf'
omni='/home/zhouchi/database/Annotation/hg38_vcf/1000G_phase1.snps.high_confidence.hg38.vcf'
OneKG='/home/zhouchi/database/Annotation/hg38_vcf/1000G_omni2.5.hg38.vcf'
mills='/home/zhouchi/database/Annotation/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf'

######somatic mutation calling#########
echo -e "\033[40;32mcalling somatic mutaiton using GATK......\t"$DATE"\033[0m"
#set -x
############################
# RNA-seq raw data mapping #
############################

###paired end###
STAR_pefq2bam(){
    local input_fq1=$1
    local input_fq2=$2
    local output_bam=$3
    local output_dup=$4
    local CPU=$5
    STAR --runThreadN $CPU --twopassMode Basic --readFilesIn $input_fq1 $input_fq2 --genomeDir $STAR_INDEX --outFileNamePrefix alignments/${PREFIX} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64424509440
    rm -rf alignments/${PREFIX}_STA* alignments/${prefix}Log* alignments/${PREFIX}SJ.out.tab
    stringtie alignments/${PREFIX}Aligned.sortedByCoord.out.bam -o expression/${PREFIX}.gtf -p 8 -G ${GTF} -A expression/${PREFIX}_gene_abund.tab
    #add_output_bam=`${iTuNES_BIN_PATH}/newtf`.bam; #mkfifo $add_output_bam
    #$AddOrReplaceReadGroups I=alignments/${PREFIX}Aligned.sortedByCoord.out.bam O=alignments/$add_output_bam SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20      # Add your own annotation
    #deduped=${PREFIX}_deduped.bam;
    #$MarkDuplicates I=alignments/$add_output_bam O=alignments/$deduped M=$output_dup REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 
    #$BuildBamIndex I=alignments/$deduped O=alignments/$deduped.bai VALIDATION_STRINGENCY=SILENT
    #$SplitNCigarReads -R $REFERENCE -I  alignments/$deduped -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o $output_bam #-fixNDN
    #rm -f alignments/$deduped alignments/$add_output_bam
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
    $BaseRecalibrator  -R $REFERENCE -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3 -o alignments/$before_table
    $PrintReads        -R $REFERENCE -I $input_bam -BQSR alignments/$before_table -o $output
    #$BaseRecalibrator  -R $REFERENCE -I $input_bam -knownSites $known_indels1 -knownSites $known_indels2 -knownSites $known_indels3  -BQSR $before_table -o $after_table
    #$AnalyzeCovariates -R $REFERENCE -before $before_table -after $after_table -plots $recal_plots -l DEBUG        
    rm alignments/$before_table
    # note that Rscript will be executed. many packages are needed. add -l DEBUG to check which library are needed 
    # if the execution failed. The following library are required:
    # reshape, gplots, ggplot2, bitops, caTools, colorspace, gtools, gData/Storage1, RColorBrewer, gsalib
}
call_variants_varscan(){
	local input_rna_bam=$1
	local out_vcf=$2
    cd somatic_mutation
    if [ -f $out_vcf ];then
    	echo "out vcf file have exist,please remove it!"
    	exit 1
   	fi
    rm -rf ${PREFIX}_tumor.fifo
	mkfifo ${PREFIX}_tumor.fifo
	samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../${input_rna_bam} > ${PREFIX}_tumor.fifo &
	java -jar /home/zhouchi/neoantigen_pipeline/varscan/VarScan.v2.4.2.jar mpileup2snp ${PREFIX}_tumor.fifo --output-vcf 1 > ../${out_vcf}
	rm ${PREFIX}_tumor.fifo
	cd ..
}
call_indel_varscan(){
	local input_rna_bam=$1
	local out_vcf=$2
    cd somatic_mutation
    rm -rf ${PREFIX}_tumor.fifo
	mkfifo ${PREFIX}_tumor.fifo
	samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../${input_rna_bam} > ${PREFIX}_tumor.fifo &
	java -jar /home/zhouchi/neoantigen_pipeline/varscan/VarScan.v2.4.2.jar mpileup2indel ${PREFIX}_tumor.fifo --output-vcf 1 > ../${out_vcf}
	rm ${PREFIX}_tumor.fifo
	cd ..	
}
strelka_indel_calling(){
	local $PREFIX
	python /usr/local/bin/configureStrelkaGermlineWorkflow.py --bam alignments/${PREFIX}_STAR_realign_recal.bam --referenceFasta=${REFERENCE} --config=/usr/local/bin/configureStrelkaGermlineWorkflow.py.ini --runDir=strelka_indel --rna
	python strelka_indel/runWorkflow.py -m local -j 8 -q mel_21_1_strelka -g 32 --quiet
}
strelkaindel2fasta(){
	local VEP_CACHE=$1
	local PREFIX=$2
	vep -i strelka_indel/results/variants/germline.indels.vcf.gz --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt -o netmhc -s ${PREFIX}_strelka
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_strelka
}
varscanindel2fasta(){
	local VEP_CACHE=$1
	local PREFIX=$2
	python ${iTuNES_BIN_PATH}/varscan_indel_preprocess.py -i somatic_mutation/${PREFIX}.indel -o somatic_mutation -s ${PREFIX}
	vep -i somatic_mutation/${PREFIX}_varscan_indel.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt -o netmhc -s ${PREFIX}_varscan
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_varscan
}
variant_calling(){
	local T_fastq_1=$1
	local T_fastq_2=$2
	local PREFIX=$3
	local CPU=$4
	STAR_pefq2bam	$T_fastq_1	$T_fastq_2	alignments/${PREFIX}_STAR.bam	alignments/${PREFIX}_STAR.dup $CPU
	#realignment	alignments/${PREFIX}_STAR.bam        $OneKG     $mills  alignments/${PREFIX}_STAR_realign.bam
	#baserecalibration       alignments/${PREFIX}_STAR_realign.bam      $dbsnp138      $OneKG     $mills   alignments/${PREFIX}_STAR_realign_recal.bam	alignments/${PREFIX}_STAR_realign_recal.pdf
	#call_variants_varscan alignments/${PREFIX}_STAR_realign_recal.bam somatic_mutation/${PREFIX}_varscan_snv.vcf
	#call_indel_varscan alignments/${PREFIX}_STAR_realign_recal.bam somatic_mutation/${PREFIX}_varscan_indel.vcf
	#strelka_indel_calling ${PREFIX}
}
genefusion(){
	local PREFIX=$1
	local CPU=$2
	cd gene_fusion
	rm -rf result
	eric_db='/home/zhouchi/database/Annotation/ericscript_db'
	eric_path=/home/zhouchi/software/ericscript-0.5.5/ericscript.pl
	$eric_path -db $eric_db --refid homo_sapiens -name ${PREFIX}_genefusion -o result $T_fastq_1 $T_fastq_2 -p $CPU
	cd ..
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
hlatyping(){
	local PREFIX=$1
	local CPU=$2
	rm -rf HLAtyping/*
	razers3 -i 95 -tc $CPU -m 1 -dr 0 -o HLAtyping/fished_test_1.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_1
	razers3 -i 95 -tc $CPU -m 1 -dr 0 -o HLAtyping/fished_test_2.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_2
	samtools bam2fq HLAtyping/fished_test_1.bam > HLAtyping/test_1_fished.fastq 
	samtools bam2fq HLAtyping/fished_test_2.bam > HLAtyping/test_2_fished.fastq
	rm HLAtyping/fished_test_1.bam HLAtyping/fished_test_2.bam
	python /home/zhouchi/software/OptiType/OptiTypePipeline.py -i HLAtyping/test_1_fished.fastq HLAtyping/test_2_fished.fastq --rna -o HLAtyping/
	rm HLAtyping/test_1_fished.fastq HLAtyping/test_2_fished.fastq
	dir=`ls HLAtyping/`
	python ${iTuNES_BIN_PATH}/optitype_ext.py -i HLAtyping/${dir}/${dir}_result.tsv -o HLAtyping -s ${PREFIX}
}

VarscanSnvVcf2Fasta2Netmhcpan2NetCTLPAN(){
	local PREFIX=$1
	local Exp_file=$2
	local Binding_Aff_Fc_Cutoff=$3
	local Binding_Aff_Cutoff=$4 
	local Fpkm_Cutoff=$5
	local VEP_CACHE=$6
	vep -i somatic_mutation/${PREFIX}_varscan_snv.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_varscan_snv_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/snv2fasta.py -i somatic_mutation/${PREFIX}_varscan_snv_vep_ann.txt -o netmhc -s ${PREFIX}_varscan
	netmhcpan netmhc/${PREFIX}_varscan_snv.fasta ${hla_str} ${PREFIX}_varscan_snv_netmhc.txt netmhc 15000
	python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_varscan_snv_netmhc.txt -g netmhc/${PREFIX}_varscan_snv.fasta -o netmhc -s ${PREFIX}_varscan_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
	python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_varscan_snv_final_neo_candidate.txt -o netctl -s ${PREFIX}_varscan_snv

}


IndelVcf2Fasta2Netmhcpan2NetCTLPAN(){
	local PREFIX=$1
	local Exp_file=$2
	local Binding_Aff_Fc_Cutoff=$3
	local Binding_Aff_Cutoff=$4 
	local Fpkm_Cutoff=$5
	local VEP_CACHE=$6

	vep -i somatic_mutation/${PREFIX}_varscan_indel.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt -o netmhc -s ${PREFIX}_varscan
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_varscan_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_varscan
	vep -i strelka_indel/results/variants/germline.indels.vcf.gz --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
	python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt -o netmhc -s ${PREFIX}_strelka
	python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i somatic_mutation/${PREFIX}_strelka_indel_vep_ann.txt  -o netmhc -s ${PREFIX}_strelka
	cat netmhc/${PREFIX}_strelka_del.fasta > netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_strelka_ins.fasta >> netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_varscan_del.fasta >> netmhc/${PREFIX}_indel.fasta
	cat netmhc/${PREFIX}_varscan_ins.fasta >> netmhc/${PREFIX}_indel.fasta
	netmhcpan netmhc/${PREFIX}_indel.fasta ${hla_str} ${PREFIX}_indel_netmhc.txt netmhc 100
	python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_indel_netmhc.txt -g netmhc/${PREFIX}_indel.fasta -o netmhc -s ${PREFIX}_indel -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
	python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_indel_final_neo_candidate.txt -o netctl -s ${PREFIX}_indel

}

GeneFusion2Fasta2Netmhcpan2NetCTLPAN(){
	local PREFIX=$1
	local GeneFusion_Cutoff=$2
	python ${iTuNES_BIN_PATH}/genefusion2fasta.py -i gene_fusion/result/${PREFIX}_genefusion.results.total.tsv -c ${GeneFusion_Cutoff} -o netmhc/ -s ${PREFIX}
	netmhcpan netmhc/${PREFIX}_gene_fusion.fasta ${hla_str} ${PREFIX}_genefusion_netmhc.txt netmhc 100
	python ${iTuNES_BIN_PATH}/gf_netMHC_result_parse.py -i netmhc/${PREFIX}_genefusion_netmhc.txt -t netmhc/${PREFIX}_gene_fusion.fasta -g gene_fusion/result/${PREFIX}_genefusion.results.total.tsv -o netmhc -s ${PREFIX} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
	python ${iTuNES_BIN_PATH}/gf_netctl.py -i netmhc/${PREFIX}_genefusion_final_neo_candidate.txt -o netctl -s ${PREFIX}_gf
}


startTime=`date +"%s.%N"`
pipeline_all(){
	set -x
	local PREFIX=$1
	echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
	if [ "$Hla_Type" = 'None' ];then
		echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
		#hlatyping $PREFIX $CPU
		#hla_str=`cat HLAtyping/${PREFIX}_optitype_hla_type`
	else
		hla_str=${Hla_Type}
	fi
	#genefusion $PREFIX $CPU 
	#GeneFusion2Fasta2Netmhcpan2NetCTLPAN $PREFIX $GeneFusion_Cutoff
	echo -e "\033[40;32mmapping......\t"$DATE"\033[0m"
	variant_calling ${T_fastq_1} ${T_fastq_2} $PREFIX $CPU
	echo -e "\033[40;32mindel calling ......\t"$DATE"\033[0m"
	echo -e "\033[40;32fasta generate ......\t"$DATE"\033[0m"
	#VEP_CACHE=/home/zhouchi/database/Annotation/vep_data
	#VarscanSnvVcf2Fasta2Netmhcpan2NetCTLPAN $PREFIX ${Exp_file} ${Binding_Aff_Fc_Cutoff} ${Binding_Aff_Cutoff} ${Fpkm_Cutoff} $VEP_CACHE
	echo -e "\033[40;32Runnning netmhcpan and netCTLPAN ......\t"$DATE"\033[0m"
	echo -e "\033[40;32Runnning pyclone ......\t"$DATE"\033[0m"
	set +x
}
pipeline_all $PREFIX


###finished
endTime=`date +"%s.%N"`
echo `awk -v x1="$(echo $endTime | cut -d '.' -f 1)" -v x2="$(echo $startTime | cut -d '.' -f 1)" -v y1="$[$(echo $endTime | cut -d '.' -f 2) / 1000]" -v y2="$[$(echo $startTime | cut -d '.' -f 2) /1000]" 'BEGIN{printf "RunTIme:%.6f s",(x1-x2)+(y1-y2)/1000000}'`
DATE=`date --date="-24 hour"`
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"
