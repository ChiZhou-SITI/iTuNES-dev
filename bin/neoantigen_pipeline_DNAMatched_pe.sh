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
TEMP=`getopt -o p:q:s:t:o:g:P:I:r:c:C:a:b:E:f: \
--long N_fq_1:,N_fq_2:,T_fq_1:,T_fq_2:,outpath:,genome:,prefix:,bwa_index:,reference:,CPU:Coverage:,binding_aff_fc_cutoff:,binding_aff_cutoff:,expression_file:,fpkm_cutoff: \
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
######check other files########



####################
###Align with bwa###
####################
###align
:<<!
DATE=`date --date="-24 hour"`
echo -e "\033[40;36;1m\033[1m---Align......---\t"$DATE"\033[0m"
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mbwa alignment......\t"$DATE"\033[0m"
set -x

bwa mem -t ${CPU} ${BWA_INDEX}/${GENOME} ${N_fastq_1} ${N_fastq_2} > alignments/tmp_${PREFIX}_normal.sam \
	2>log_file/${PREFIX}_normal_bwa_aln.log 
	
bwa mem -t ${CPU} ${BWA_INDEX}/${GENOME} ${T_fastq_1} ${T_fastq_2} > alignments/tmp_${PREFIX}_tumor.sam \
	2>log_file/${PREFIX}_tumor_bwa_aln.log

samtools view -bhS alignments/tmp_${PREFIX}_normal.sam > alignments/tmp_${PREFIX}_normal.bam 
samtools view -bhS alignments/tmp_${PREFIX}_tumor.sam > alignments/tmp_${PREFIX}_tumor.bam

samtools sort -@ ${CPU} -m 8G alignments/tmp_${PREFIX}_normal.bam alignments/${PREFIX}_normal_unfilter 
samtools sort -@ ${CPU} -m 8G alignments/tmp_${PREFIX}_tumor.bam alignments/${PREFIX}_tumor_unfilter

rm alignments/tmp_${PREFIX}_normal.bam alignments/tmp_${PREFIX}_tumor.bam alignments/tmp_${PREFIX}_normal.sam alignments/tmp_${PREFIX}_tumor.sam
samtools index alignments/${PREFIX}_normal_unfilter.bam 
samtools index alignments/${PREFIX}_tumor_unfilter.bam 

samtools view -b alignments/${PREFIX}_normal_unfilter.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > alignments/${PREFIX}_normal.bam 
samtools view -b alignments/${PREFIX}_tumor_unfilter.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > alignments/${PREFIX}_tumor.bam
samtools index alignments/${PREFIX}_normal.bam 
samtools index alignments/${PREFIX}_tumor.bam
rm alignments/${PREFIX}_normal_unfilter.bam alignments/${PREFIX}_tumor_unfilter.bam alignments/${PREFIX}_normal_unfilter.bam.bai alignments/${PREFIX}_tumor_unfilter.bam.bai
set +x
!
#######################
###Filter alignments###
#######################
###mark duplicates with picard
DATE=`date --date="-24 hour"`
:<<!
echo -e "\033[40;32mrun picard mark duplicates on non-UMI......\t"$DATE"\033[0m"
set -x
time java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar MarkDuplicates \
	INPUT=alignments/${PREFIX}_normal.bam OUTPUT=alignments/${PREFIX}_normal_marked.bam \
	METRICS_FILE=${PREFIX}_normal_dup_qc.txt ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
#	READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
time java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar MarkDuplicates \
	INPUT=alignments/${PREFIX}_tumor.bam OUTPUT=alignments/${PREFIX}_tumor_marked.bam \
	METRICS_FILE=${PREFIX}_tumor_dup_qc.txt ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT
set +x
mv ${PREFIX}_normal_dup_qc.txt bamstat
mv ${PREFIX}_tumor_dup_qc.txt bamstat
###fliter bam on flags and threashold
DATE=`date --date="-24 hour"`
echo -e "\033[40;32mfilter on flag 512 and threadshold score: "$SCORE"......\t"$DATE"\033[0m"
samtools view -F 512 -b alignments/${PREFIX}_normal_marked.bam > alignments/${PREFIX}_normal_mkdup_filter.bam
samtools view -F 512 -b alignments/${PREFIX}_tumor_marked.bam > alignments/${PREFIX}_tumor_mkdup_filter.bam
rm alignments/${PREFIX}_normal_marked.bam alignments/${PREFIX}_tumor_marked.bam
!
:<<!
######somatic mutation calling#########
echo -e "\033[40;32mcalling somatic mutaiton using varscan2......\t"$DATE"\033[0m"
set -x
cd somatic_mutation
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../alignments/${PREFIX}_normal_mkdup_filter.bam > ${PREFIX}_normal.fifo &
samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../alignments/${PREFIX}_tumor_mkdup_filter.bam > ${PREFIX}_tumor.fifo &
java -jar /home/zhouchi/neoantigen_pipeline/varscan/VarScan.v2.4.2.jar somatic ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX}
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
java -jar /home/zhouchi/neoantigen_pipeline/varscan/VarScan.v2.4.2.jar processSomatic ${PREFIX}.snp
cd ..
set +x
!
##########neoantigen identification#####
#######HLA TYPING#######
echo -e "\033[40;32mHLAtyping......\t"$DATE"\033[0m"
#######optitype#############
:<<!
set -x
razers3 -i 95 -tc 8 -m 1 -dr 0 -o HLAtyping/fished_test_1.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_1
razers3 -i 95 -tc 8 -m 1 -dr 0 -o HLAtyping/fished_test_2.bam /home/zhouchi/software/OptiType/data/hla_reference_dna.fasta $T_fastq_2
samtools bam2fq HLAtyping/fished_test_1.bam > HLAtyping/test_1_fished.fastq 
samtools bam2fq HLAtyping/fished_test_2.bam > HLAtyping/test_2_fished.fastq
rm HLAtyping/fished_test_1.bam HLAtyping/fished_test_2.bam
python /home/zhouchi/software/OptiType/OptiTypePipeline.py -i $T_fastq_1 $T_fastq_2 --dna -o HLAtyping/
rm HLAtyping/test_1_fished.fastq HLAtyping/test_2_fished.fastq
if [ -f HLAtyping/${PREFIX}_optitype_hla_type ];then
	rm HLAtyping/${PREFIX}_optitype_hla_type
fi
dir=`ls HLAtyping/`
python ${iTuNES_BIN_PATH}/optitype_ext.py -i HLAtyping/${dir}/${dir}_result.tsv -o HLAtyping -s ${PREFIX}
rm -rf HLAtyping/${dir}
############################
set +x
!
set -x
#######SNV derived neoantigens########
echo -e "\033[40;32mVEP Annotation......\t"$DATE"\033[0m"
VEP_CACHE=/home/zhouchi/database/Annotation/vep_data
sed '1d' somatic_mutation/${PREFIX}.snp.Somatic | awk -F '\t' '{print $1,$2,$2,$3"/"$4}' > somatic_mutation/${PREFIX}_vep_input.vcf
VEP_CACHE=/data/PUBLIC/VEP_DATA
vep -i somatic_mutation/${PREFIX}_vep_input.vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o somatic_mutation/${PREFIX}_snv_vep_ann.txt --force_overwrite
set +x
set -x
###indel identification###
echo -e "\033[40;32mindel calling ......\t"$DATE"\033[0m"
:<<! 
samtools index alignments/${PREFIX}_normal_mkdup_filter.bam
samtools index alignments/${PREFIX}_tumor_mkdup_filter.bam
!
cd somatic_mutation
#DEL_PREFIX=${PREFIX}_del
#INS_PREFIX=${PREFIX}_ins
#DATE=`date --date="-24 hour"`
:<<!
echo -e "\033[40;32mfinish indel calling and annotation ......\t"$DATE"\033[0m"
echo -e ../alignments/${PREFIX}_tumor_mkdup_filter.bam'\t250\t'${PREFIX} > ${PREFIX}_config.txt
echo -e ../alignments/${PREFIX}_normal_mkdup_filter.bam'\t250\t'${PREFIX} >> ${PREFIX}_config.txt
CHR_ARR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
for CHR in ${CHR_ARR[*]}
do
{
	pindel -f /home/zhouchi/database/Annotation/Fasta/human.fasta -i ${PREFIX}_config.txt -c $CHR -o ${PREFIX}_${CHR} -E 0.90 -T 32
}
done
DATE1=`date --date="-24 hour"`
echo -e "\033[40;32mfinish......\t"$DATE1"\033[0m"
for CHR in ${CHR_ARR[*]}
do
{
	cat ${PREFIX}_${CHR}_D >> ${PREFIX}_D
	cat ${PREFIX}_${CHR}_SI >> ${PREFIX}_SI
	rm ${PREFIX}_${CHR}_D
	rm ${PREFIX}_${CHR}_SI
}
done
!
VEP_CACHE=/home/zhouchi/database/Annotation/vep_data
pindel2vcf -p ${PREFIX}_D -r /home/zhouchi/database/Annotation/Fasta/human.fasta -R ucsc_hg19 -d 20090201 -v ${PREFIX}_D.vcf
vep -i ${PREFIX}_D.vcf --format vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is inframe_deletion" -o ${PREFIX}_D_id_vep_ann.txt --force_overwrite
vep -i ${PREFIX}_D.vcf --format vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is frameshift_variant" -o ${PREFIX}_D_fs_vep_ann.txt --force_overwrite
pindel2vcf -p ${PREFIX}_SI -r /home/zhouchi/database/Annotation/Fasta/human.fasta -R ucsc_hg19 -d 20090201 -v ${PREFIX}_SI.vcf
vep -i ${PREFIX}_SI.vcf --format vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is frameshift_variant" -o ${PREFIX}_SI_fs_vep_ann.txt --force_overwrite
vep -i ${PREFIX}_SI.vcf --format vcf --cache --dir $VEP_CACHE --dir_cache $VEP_CACHE --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is inframe_insertion" -o ${PREFIX}_SI_ii_vep_ann.txt --force_overwrite
cd ..
set +x
set -x
########generate fasta file########
python ${iTuNES_BIN_PATH}/snv2fasta.py -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -o netmhc -s ${PREFIX}
python ${iTuNES_BIN_PATH}/deletion2fasta.py -i somatic_mutation/${PREFIX}_D_id_vep_ann.txt -c somatic_mutation/${PREFIX}_D.vcf -o netmhc -s ${PREFIX}_id
python ${iTuNES_BIN_PATH}/deletion2fasta.py -i somatic_mutation/${PREFIX}_D_fs_vep_ann.txt -c somatic_mutation/${PREFIX}_D.vcf -o netmhc -s ${PREFIX}_fs
python ${iTuNES_BIN_PATH}/insertion2fasta.py -i somatic_mutation/${PREFIX}_SI_ii_vep_ann.txt -c somatic_mutation/${PREFIX}_SI.vcf -o netmhc -s ${PREFIX}_ii
python ${iTuNES_BIN_PATH}/insertion2fasta.py -i somatic_mutation/${PREFIX}_SI_fs_vep_ann.txt -c somatic_mutation/${PREFIX}_SI.vcf -o netmhc -s ${PREFIX}_fs
########run netMHCpan###########
#netMHCpan -a $hla_str -f netmhc/${PREFIX}_snv.fasta -l 8,9,10,11 > netmhc/${PREFIX}_snv_netmhc.txt
#####snv######
netmhcpan(){
	local input_fasta=$1
	local hla=$2
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
	for file in $filelist
	do
	{
		netMHCpan -a ${hla} -f ${out_dir}/tmp/${file} -l 8,9,10,11 > ${out_dir}/tmp/${file}_tmp_netmhc.txt
		rm ${out_dir}/tmp/${file}
	} &
	done
	wait
	filelist1=`ls ${out_dir}/tmp/`
	for file_r in $filelist1
	do
	{
		cat ${out_dir}/tmp/$file_r >> ${out_dir}/${netmhc_out}
		rm ${out_dir}/tmp/${file_r}	
	}
	done
	rm -rf 	${out_dir}/tmp
}
hla_str=`cat HLAtyping/${PREFIX}_optitype_hla_type`
netmhcpan netmhc/${PREFIX}_snv.fasta ${hla_str} ${PREFIX}_snv_netmhc.txt netmhc 1000
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_snv_netmhc.txt -g netmhc/${PREFIX}_snv.fasta -o netmhc -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
netmhcpan netmhc/${PREFIX}_id_del.fasta ${hla_str} ${PREFIX}_id_del_netmhc.txt netmhc 6000
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_id_del_netmhc.txt -g netmhc/${PREFIX}_id_del.fasta -o netmhc -s ${PREFIX}_id_del -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
netmhcpan netmhc/${PREFIX}_fs_del.fasta ${hla_str} ${PREFIX}_fs_del_netmhc.txt netmhc 50000
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_fs_del_netmhc.txt -g netmhc/${PREFIX}_fs_del.fasta -o netmhc -s ${PREFIX}_fs_del -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
netmhcpan netmhc/${PREFIX}_ii_ins.fasta ${hla_str} ${PREFIX}_ii_ins_netmhc.txt netmhc 8000
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_ii_ins_netmhc.txt -g netmhc/${PREFIX}_ii_ins.fasta -o netmhc -s ${PREFIX}_ii_ins -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
netmhcpan netmhc/${PREFIX}_fs_ins.fasta ${hla_str} ${PREFIX}_fs_ins_netmhc.txt netmhc 8000
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i netmhc/${PREFIX}_fs_ins_netmhc.txt -g netmhc/${PREFIX}_fs_ins.fasta -o netmhc -s ${PREFIX}_fs_ins -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff}
set +x
#######netCTLpan########
set -x
echo -e "\033[40;32mnetCTLpan prediction\t"$DATE"\033[0m"
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_snv_final_neo_candidate.txt -o netctl -s ${PREFIX}_snv
if [ "$?" -ne 0 ]; then echo "Running Snv netCTLPAN failed!"; exit 1; else echo "Running Snv netCTLPAN finished!"; fi 
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_id_del_final_neo_candidate.txt -o netctl -s ${PREFIX}_id_del
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_fs_del_final_neo_candidate.txt -o netctl -s ${PREFIX}_fs_del
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_ii_ins_final_neo_candidate.txt -o netctl -s ${PREFIX}_ii_ins
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i netmhc/${PREFIX}_fs_ins_final_neo_candidate.txt -o netctl -s ${PREFIX}_fs_ins
echo -e "\033[40;32mFinish netCTLpan\t"$DATE"\033[0m"
set +x
echo -e "\033[40;32mGene fusion detection\t"$DATE"\033[0m"
####################################################
########## detect gene fusion using DNA-seq#########
####################################################
######detect copynumber profile#########
echo -e "\033[40;32mdetect copy number profile......\t"$DATE"\033[0m"
:<<!
set -x
cd copynumber_profile
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../alignments/${PREFIX}_normal_mkdup_filter.bam > ${PREFIX}_normal.fifo &
samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ../alignments/${PREFIX}_tumor_mkdup_filter.bam > ${PREFIX}_tumor.fifo &
java -jar /home/zhouchi/neoantigen_pipeline/varscan/VarScan.v2.4.2.jar copynumber ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX}
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
set +x
!
set -x
Rscript ${iTuNES_BIN_PATH}/sequenza_test.R somatic_mutation/${PREFIX}.snp copynumber_profile/${PREFIX}.copynumber copynumber_profile/ ${PREFIX}
if [ "$?" -ne 0 ]; then echo "snv netMHCMHC result parsing failed"; exit 1; else echo "snv netMHCMHC result parsing finished"; fi 
#echo $COVERAGE
python ${iTuNES_BIN_PATH}/pyclone_input.py -n netctl/${PREFIX}_snv_netctl_concact.txt -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -s somatic_mutation/${PREFIX}.snp.Somatic -c copynumber_profile/${PREFIX}_seg_copynumber.txt -o copynumber_profile -S ${PREFIX} -C ${COVERAGE}
if [ "$?" -ne 0 ]; then echo "prepare pyclone input failed"; exit 1; else echo "Prepare pyclone input finished!"; fi 
TUMOR_CONTENT=`cat copynumber_profile/${PREFIX}_cellularity.txt`
PyClone run_analysis_pipeline --in_files copynumber_profile/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir pyclone 
if [ "$?" -ne 0 ]; then echo "Pyclone running failed"; exit 1; else echo "Pyclone running finished!"; fi 
python ${iTuNES_BIN_PATH}/neo_pyclone_annotation.py -n netctl/${PREFIX}_snv_netctl_concact.txt -i somatic_mutation/${PREFIX}_snv_vep_ann.txt -s pyclone/tables/loci.tsv -o netctl -S ${PREFIX}
if [ "$?" -ne 0 ]; then echo "Pyclone annotation failed"; exit 1; else echo "Pyclone annotation finished!"; fi 
set +x

######
###finished
echo -e "\033[40;32;1mfinished\t"$DATE"\033[0m"