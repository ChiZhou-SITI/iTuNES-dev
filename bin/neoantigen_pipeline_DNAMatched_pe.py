import os,sys,time
import multiprocessing
import threading
import Queue
import shutil 
import Queue
q=Queue.Queue()
import subprocess
import getopt

global CPU
global optitype_ref_path
global opitype_fold
global opitype_ext
global BWA_INDEX
global GENOME
global prefix
global REFERENCE
global opitype_out_fold


prefix='mel_21_1'
CPU=8
GENOME='human'
BWA_INDEX='/home/zhouchi/database/Annotation/Index/bwa/'+GENOME
REFERENCE='/home/zhouchi/database/Annotation/Fasta/'+GENOME+'.fasta'
VEP_CACHE='/home/zhouchi/database/Annotation/vep_data/'
opitype_fold='/home/zhouchi/software/OptiType/OptiTypePipeline.py'
tumor_fastq_path_first='/home/zhouchi/raw_data/phs001005/fastq/SRR2672887/SRR2672887_1.fastq'
tumor_fastq_path_second='/home/zhouchi/raw_data/phs001005/fastq/SRR2672887/SRR2672887_2.fastq'
normal_fastq_path_first='/home/zhouchi/raw_data/phs001005/fastq/SRR2672888/SRR2672888_1.fastq'
normal_fastq_path_second='/home/zhouchi/raw_data/phs001005/fastq/SRR2672888/SRR2672888_2.fastq'
opitype_out_fold='HLAtyping'
alignment_out_fold='alignments'
netmhc_out_fold='netmhc'
somatic_mutation_fold='somatic_mutation'
strelka_out_fold='strelka_indel'
pyclone_fold='pyclone'
snv_fasta_file=prefix+'_snv.fasta'
hla_str='HLA-A02:01'
snv_netmhc_out_file=prefix+'_snv_netmhc.txt'
indel_fasta_file=prefix+'_indel.fasta'
indel_netmhc_out_file=prefix+'_indel_netmhc.txt'
split_num=1000
exp_file='no_exp'
binding_fc_aff_cutoff=1
binding_aff_cutoff=500
fpkm_cutoff=1
coverage=1
netctl_out_fold='netctl'
if not os.path.exists(netctl_out_fold):
	os.mkdir(netctl_out_fold)
opitype_ext='${iTuNES_BIN_PATH}/optitype_ext.py'
java_picard_markdup='java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar MarkDuplicates'
GATK="/home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar"
AddOrReplaceReadGroups='java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups'
BuildBamIndex="java -Xmx4G -jar /home/zhouchi/software/gatk_pre/picard-tools-2.3.0/picard.jar BuildBamIndex"
RealignerTargetCreator="java -Xmx4G -jar /home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 8 -dt NONE"
IndelRealigner="java -Xmx4G -jar /home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar -T IndelRealigner"
BaseRecalibrator="java -Xmx4G -jar /home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8"
PrintReads="java -Xmx4G -jar /home/zhouchi/software/gatk_pre/GenomeAnalysisTK.jar -T PrintReads -nct 8 -dt NONE"
dbsnp138='/home/zhouchi/database/Annotation/hg38_vcf/dbsnp_138.hg38.vcf'
hapmap='/home/zhouchi/database/Annotation/hg38_vcf/hapmap_3.3.hg38.vcf'
omni='/home/zhouchi/database/Annotation/hg38_vcf/1000G_phase1.snps.high_confidence.hg38.vcf'
OneKG='/home/zhouchi/database/Annotation/hg38_vcf/1000G_omni2.5.hg38.vcf'
mills='/home/zhouchi/database/Annotation/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf'
cosmic='/home/zhouchi/database/Annotation/hg38_vcf/CosmicCodingMuts_chr_M_sorted.vcf'
if not os.path.exists(alignment_out_fold):
	os.mkdir(alignment_out_fold)
else:
	shutil.rmtree(alignment_out_fold)  
	os.mkdir(alignment_out_fold) 
	print 'fold already exist'
def hlatyping(raw_fastq_path_first,raw_fastq_path_second):
	if not os.path.exists(opitype_out_fold):
		os.mkdir(opitype_out_fold)
	else:
		shutil.rmtree(opitype_out_fold)  
		os.mkdir(opitype_out_fold) 
		print 'fold already exist'
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' --dna -o ' + opitype_out_fold
	print cmd_hla
	os.system(cmd_hla)
	result_dir=os.listdir(opitype_out_fold)
	print result_dir[0]
	hla_result_path=opitype_out_fold+'/'+result_dir[0]+'/'+result_dir[0]+'_result.tsv'
	print hla_result_path
	cmd_hla_ext = 'python ' + opitype_ext + ' -i ' + hla_result_path + ' -o ' + opitype_out_fold + ' -s ' + prefix
	print cmd_hla_ext
	os.system(cmd_hla_ext)
	print 'hla type process done.'

def mapping_qc_gatk_preprocess(fastq_1_path,fastq_2_path,fastq_type):
	cmd_bwa='bwa mem -t '+ str(CPU) + ' ' + BWA_INDEX + '/' + GENOME + ' ' + fastq_1_path + ' ' +fastq_2_path + ' > ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam'
	cmd_samtools_1='samtools view -bhS -@ '+ str(CPU) + ' ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam' + ' > ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam'
	cmd_samtools_sort='samtools sort -@ ' + '4' + ' -m 8G ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam' + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter'
	cmd_samtools_index_1='samtools index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter.bam'
	cmd_select_chr='samtools view -b -@ ' + str(CPU) + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter.bam' + ' chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ' + alignment_out_fold+'/' + prefix + '_'+fastq_type+'.bam'
	cmd_samtools_2='samtools index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam'
	cmd_picard=java_picard_markdup + ' INPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_'+fastq_type+'_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT'
	cmd_samtools_index_2='samtools index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam'
	cmd_add_readgroup=AddOrReplaceReadGroups + ' I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20'
	cmd_buildbamindex=BuildBamIndex + ' I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam.bai' + ' VALIDATION_STRINGENCY=SILENT'
	cmd_RealignerTargetCreator=RealignerTargetCreator + ' -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.intervals'
	cmd_IndelRealigner=IndelRealigner + ' -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -targetIntervals ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.intervals' + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam'
	cmd_BaseRecalibrator=BaseRecalibrator + ' -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table'
	cmd_PrintReads=PrintReads + ' -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '_recal.bam'
	#print cmd_bwa
	#print cmd_samtools_1
	#print cmd_samtools_sort
	#print cmd_samtools_index_1
	#print cmd_select_chr
	#print cmd_samtools_2
	#print cmd_picard
	#print cmd_samtools_index_2
	#print cmd_add_readgroup
	#print cmd_buildbamindex
	#print cmd_RealignerTargetCreator
	#print cmd_IndelRealigner
	#print cmd_BaseRecalibrator
	#print cmd_PrintReads
	os.system(cmd_bwa)
	os.system(cmd_samtools_1)
	os.system(cmd_samtools_sort)
	os.system(cmd_samtools_index_1)
	os.system(cmd_select_chr)
	os.system(cmd_samtools_2)
	os.system(cmd_picard)
	os.system(cmd_samtools_index_2)
	os.system(cmd_add_readgroup)
	os.system(cmd_buildbamindex)
	os.system(cmd_RealignerTargetCreator)
	os.system(cmd_IndelRealigner)
	os.system(cmd_BaseRecalibrator)
	os.system(cmd_PrintReads)
def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num):
	str_proc=r'''
set -x
input_fasta=%s
hla_str=%s
netmhc_out=%s
out_dir=%s
split_num=%s
if [ -d ${out_dir}/tmp ];then
	rm -rf ${out_dir}/tmp
	mkdir -p ${out_dir}/tmp
else
	mkdir -p ${out_dir}/tmp
fi
if [ -f ${out_dir}/${netmhc_out} ];then
	rm ${out_dir}/${netmhc_out}
fi
split -l ${split_num} ${out_dir}/${input_fasta} ${out_dir}/tmp/
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
set +x
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_somatic_caling_drift(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,snv_fasta_file,hla_str,snv_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold):
	str_proc = r'''
set -e
somat_f=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
if [ ! -d ${somat_f} ];then
	mkdir ${somat_f}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
rm -rf ${somat_f}/*
cd ${somat_f}
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ../${alignment_fold}/${PREFIX}_normal_recal.bam  > ${PREFIX}_normal.fifo &
samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ../${alignment_fold}/${PREFIX}_tumor_recal.bam  > ${PREFIX}_tumor.fifo &
java -jar /home/zhouchi/software/varscan/VarScan.v2.4.2.jar somatic ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX} #--output-vcf 1
java -jar /home/zhouchi/software/varscan/VarScan.v2.4.2.jar processSomatic ${PREFIX}.snp
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
sed '1d' ${somat_f}/${PREFIX}.snp.Somatic | awk -F '\t' '{print $1,$2,$2,$3"/"$4}' > ${somat_f}/${PREFIX}_snv_vep_input.vcf
vep -i ${somat_f}/${PREFIX}_snv_vep_input.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol --offline -o ${somat_f}/${PREFIX}_snv_all_vep_ann.txt
vep -i ${somat_f}/${PREFIX}_snv_vep_input.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o ${somat_f}/${PREFIX}_snv_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/snv2fasta.py -i ${somat_f}/${PREFIX}_snv_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
'''%(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
	netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num)
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.txt -g ${netmhc_out}/${PREFIX}_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_snv
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')

def indel_calling_drift(somatic_mutation_fold,strelka_out_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold):
	str_proc1='''
somatic_mutation=%s
PREFIX=%s
vep_cache=%s
netmhc_out=%s
python ${iTuNES_BIN_PATH}/varscan_indel_preprocess.py -i ${somatic_mutation}/${PREFIX}.indel -o ${somatic_mutation} -s ${PREFIX}
vep -i ${somatic_mutation}/${PREFIX}_varscan_indel.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(somatic_mutation_fold,PREFIX,vep_cache,netmhc_out_fold)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')
	str_proc2=r'''
strelka_fold=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
cpu=%s
if [ -d ${strelka_fold} ];then
	rm -rf ${strelka_fold}
fi
python /usr/local/bin/configureStrelkaSomaticWorkflow.py --tumorBam=${alignment_fold}/${PREFIX}_tumor_recal.bam --normalBam=${alignment_fold}/${PREFIX}_normal_recal.bam --referenceFasta=${REFERENCE} --config=/usr/local/bin/configureStrelkaSomaticWorkflow.py.ini --runDir=${strelka_fold} --exome
python ${strelka_fold}/runWorkflow.py -m local -j $cpu -q ${PREFIX}_strelka -g 32 --quiet
vep -i ${strelka_fold}/results/variants/somatic.indels.vcf.gz --cache --dir ${vep_cache} --dir_cache ${vep_cache} --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_strelka
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_strelka
'''%(strelka_out_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU)
	print str_proc2
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')
	str_proc3=r'''
PREFIX=%s
netmhc_out=%s
cat ${netmhc_out}/${PREFIX}_strelka_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_strelka_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_del.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
'''%(prefix,netmhc_out_fold)
	subprocess.call(str_proc3, shell=True, executable='/bin/bash')
	netMHCpan(indel_fasta_file,hla_str,indel_netmhc_out_file,netmhc_out_fold,split_num)
	str_proc4=r'''
set -e
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_indel_netmhc.txt -g ${netmhc_out}/${PREFIX}_indel.fasta -o ${netmhc_out} -s ${PREFIX}_indel -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_indel_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_indel
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	#print str_proc4
	subprocess.call(str_proc4, shell=True, executable='/bin/bash')
	
	
def pyclone_copynumber_preprocess(somatic_mutation_fold,prefix,alignment_out_fold,REFERENCE,pyclone_fold):
	str_proc=r'''
set -x
somatic_mutation=%s
PREFIX=%s
alignments=%s
REFERENCE=%s
pyclone=%s
if [ ! -d ${pyclone} ];then
	mkdir ${pyclone}
fi
chrlist=(`sed '1d' ${somatic_mutation}/${PREFIX}.snp | cut -f1 | uniq`)
chrlist_1=${chrlist[@]:0:6}
chrlist_2=${chrlist[@]:6:6}
chrlist_3=${chrlist[@]:12:6}
chrlist_4=${chrlist[@]:18:${#chrlist[@]}}
echo 'chromosome	position	base.ref	depth.normal	depth.tumor	depth.ratio	Af	Bf	zygosity.normal	GC.percent	good.reads	AB.normal	AB.tumor	tumor.strand' > ${pyclone}/${PREFIX}.seqz
for s in $chrlist_1
do
{
	str_chr='-C '$s
	python ${iTuNES_BIN_PATH}/sequenza-utils.py bam2seqz -gc /home/zhouchi/database/Annotation/human.gc50Base.txt.gz $str_chr --fasta $REFERENCE -n ${alignments}/${PREFIX}_normal_recal.bam -t ${alignments}/${PREFIX}_tumor_recal.bam | sed '1d' >> ${pyclone}/${s}.seqz
} &
done
wait
for s in $chrlist_2
do
{
	str_chr='-C '$s
	python ${iTuNES_BIN_PATH}/sequenza-utils.py bam2seqz -gc /home/zhouchi/database/Annotation/human.gc50Base.txt.gz $str_chr --fasta $REFERENCE -n ${alignments}/${PREFIX}_normal_recal.bam -t ${alignments}/${PREFIX}_tumor_recal.bam | sed '1d' >> ${pyclone}/${s}.seqz
} &
done
wait
for s in $chrlist_3
do
{
	str_chr='-C '$s
	python ${iTuNES_BIN_PATH}/sequenza-utils.py bam2seqz -gc /home/zhouchi/database/Annotation/human.gc50Base.txt.gz $str_chr --fasta $REFERENCE -n ${alignments}/${PREFIX}_normal_recal.bam -t ${alignments}/${PREFIX}_tumor_recal.bam | sed '1d' >> ${pyclone}/${s}.seqz
} &
done
wait
for s in $chrlist_4
do
{
	str_chr='-C '$s
	python ${iTuNES_BIN_PATH}/sequenza-utils.py bam2seqz -gc /home/zhouchi/database/Annotation/human.gc50Base.txt.gz $str_chr --fasta $REFERENCE -n ${alignments}/${PREFIX}_normal_recal.bam -t ${alignments}/${PREFIX}_tumor_recal.bam | sed '1d' >> ${pyclone}/${s}.seqz
} &
done
wait

for s in $chrlist
do
{
	cat ${pyclone}/${s}.seqz >> ${pyclone}/${PREFIX}.seqz
	rm ${pyclone}/${s}.seqz
}
done

gzip -f ${pyclone}/${PREFIX}.seqz > ${pyclone}/${PREFIX}.seqz.gz
python ${iTuNES_BIN_PATH}/sequenza-utils.py seqz-binning -w 50 -s ${pyclone}/${PREFIX}.seqz.gz | gzip > ${pyclone}/${PREFIX}.small.seqz.gz
Rscript ${iTuNES_BIN_PATH}/sequenza_test.R ${pyclone}/${PREFIX}.small.seqz.gz ${pyclone}/ ${PREFIX}
set +x
'''%(somatic_mutation_fold,prefix,alignment_out_fold,REFERENCE,pyclone_fold)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def pyclone_annotation(somatic_mutation_fold,prefix,pyclone_fold,netctl_fold,coverage):
	str_proc=r'''
somatic_mutation=%s
PREFIX=%s
pyclone=%s
netctl=%s
COVERAGE=%d
python ${iTuNES_BIN_PATH}/pyclone_input.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann.txt -s ${somatic_mutation}/${PREFIX}.snp.Somatic -c ${pyclone}/${PREFIX}_seg_copynumber.txt -o ${pyclone} -S ${PREFIX} -C ${COVERAGE}
TUMOR_CONTENT=`cat ${pyclone}/${PREFIX}_cellularity.txt`
PyClone run_analysis_pipeline --in_files ${pyclone}/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir ${pyclone}
python ${iTuNES_BIN_PATH}/neo_pyclone_annotation.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann.txt -s ${pyclone}/tables/loci.tsv -o ${netctl} -S ${PREFIX}
'''%(somatic_mutation_fold,prefix,pyclone_fold,netctl_fold,coverage)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

if __name__ == '__main__':
	print "start stage 1"
	#processes_1=[]
	#d1=multiprocessing.Process(target=hlatyping,args=(tumor_fastq_path_first,tumor_fastq_path_second,))
 	#processes_1.append(d1)
 	#q.put('hlatyping')
 	#d2=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(normal_fastq_path_first,normal_fastq_path_second,'normal',))
 	#processes_1.append(d2)
 	#q.put('normal_qc')
 	#d3=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(tumor_fastq_path_first,tumor_fastq_path_second,'tumor',))
 	#processes_1.append(d3)
 	#q.put('tumor_qc')
 	#for p in processes_1:
	#	p.daemon = True
	#	p.start()
	#for p in processes_1:
	#	p.join()
	print 'stage 1 done.'
	print 'start stage 2'
	#processes_2=[]
	#h1=multiprocessing.Process(target=varscan_somatic_caling_drift,args=(somatic_mutation_fold,alignment_out_fold,prefix,REFERENCE,VEP_CACHE,netmhc_out_fold,snv_fasta_file,hla_str,snv_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold))
	#processes_2.append(h1)
	#h2=multiprocessing.Process(target=indel_calling_drift,args=(somatic_mutation_fold,strelka_out_fold,alignment_out_fold,prefix,REFERENCE,VEP_CACHE,netmhc_out_fold,CPU,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold))
	#processes_2.append(h2)
	#h3=multiprocessing.Process(target=pyclone_copynumber_preprocess,args=(somatic_mutation_fold,prefix,alignment_out_fold,REFERENCE,pyclone_fold,))
	#processes_2.append(h3)	
	#for p in processes_2:
	#	p.daemon = True
	#	p.start()
	#for p in processes_2:
	#	p.join()
	print 'start stage 4'
	#processes_4=[]
	#l1=multiprocessing.Process(target=pyclone_annotation,args=(somatic_mutation_fold,prefix,pyclone_fold,netctl_out_fold,coverage))
	#processes_4.append(l1)	
	#for p in processes_4:
	#	p.daemon = True
	#	p.start()
	#for p in processes_4:
	#	p.join()













