import os,sys,time
import multiprocessing
import shutil 
import subprocess
import pandas as pd
import math
from pyper import *

def read_trimmomatic(raw_fastq_path_first,raw_fastq_path_second,trimmomatic_path,adapter_path,fastq_clean_first,fastq_clean_second,fastq_unpaired_first,fastq_unpaired_second,logfile_fold):
	cmd_trimmomatic="java -jar " + trimmomatic_path + " PE -phred33 " + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' ' + fastq_clean_first + ' ' + fastq_unpaired_first + ' ' + fastq_clean_second +  ' ' + fastq_unpaired_second + " ILLUMINACLIP:" + adapter_path + ':2:30:10' + ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 > ' + logfile_fold + '/' + fastq_type + '_trimmomatic.log' + ' 2>&1'
	print cmd_trimmomatic
	os.system(cmd_trimmomatic)

def hlatyping(raw_fastq_path_first,raw_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix):
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

def mapping_qc_gatk_preprocess(fastq_1_path,fastq_2_path,fastq_type,CPU,BWA_INDEX,GENOME,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138,OneKG,mills,logfile_fold):
	cmd_bwa=bwa_path + ' mem -t '+ str(CPU) + ' ' + BWA_INDEX + '/' + GENOME + ' ' + fastq_1_path + ' ' +fastq_2_path + ' > ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam > ' + logfile_fold + '/' + fastq_type + '_bwa_aln.log' + ' 2>&1'
	cmd_samtools_1=samtools_path + ' view -bhS -@ '+ str(CPU) + ' ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.sam' + ' > ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam'
	cmd_samtools_sort=samtools_path + ' sort -@ ' + str(CPU) + ' -m 2G ' + alignment_out_fold+'/'+'tmp_'+ prefix +'_'+fastq_type+'.bam' + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter'
	cmd_samtools_index_1=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter.bam'
	cmd_select_chr=samtools_path + ' view -b -@ ' + str(CPU) + ' ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_unfilter.bam' + ' chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ' + alignment_out_fold+'/' + prefix + '_'+fastq_type+'.bam'
	cmd_samtools_2=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam'
	cmd_picard="java -Xmx4G -jar " + java_picard_path + ' MarkDuplicates INPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_'+fastq_type+'_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT > ' + logfile_fold + '/' + fastq_type + '_markdup.log' + ' 2>&1'
	cmd_samtools_index_2=samtools_path + ' index ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam'
	cmd_add_readgroup="java -Xmx4G -jar " + java_picard_path + ' AddOrReplaceReadGroups I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20 > ' + logfile_fold + '/' + fastq_type + '_addreadgroup.log' + ' 2>&1'
	cmd_buildbamindex="java -Xmx4G -jar " + java_picard_path + ' BuildBamIndex I=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam.bai' + ' VALIDATION_STRINGENCY=SILENT > ' + buildindex_logfile_path + ' 2>&1'
	cmd_RealignerTargetCreator="java -Xmx4G -jar " + GATK_path + ' -T RealignerTargetCreator -nt 8 -dt NONE -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.intervals > ' + logfile_fold + '/' + fastq_type + '_buildbamindex.log' + ' 2>&1'
	cmd_IndelRealigner="java -Xmx4G -jar " + GATK_path + ' -T IndelRealigner -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -targetIntervals ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'.intervals' + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam > ' + logfile_fold + '/' + fastq_type + '_IndelRealigner.log' + ' 2>&1'
	cmd_BaseRecalibrator="java -Xmx4G -jar " + GATK_path + ' -T BaseRecalibrator -nct 8 -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table > ' + logfile_fold + '/' + fastq_type + '_BaseRecalibrator.log' + ' 2>&1'
	cmd_PrintReads="java -Xmx4G -jar " + GATK_path + ' -T PrintReads -nct 8 -dt NONE -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_'+fastq_type+'_mkdup_add_realign.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_'+fastq_type + '_recal.bam > ' + logfile_fold + '/' + fastq_type + '_PrintRead.log' + ' 2>&1'
	print cmd_bwa
	os.system(cmd_bwa)
	print cmd_samtools_1
	os.system(cmd_samtools_1)
	print cmd_samtools_sort
	os.system(cmd_samtools_sort)
	print cmd_samtools_index_1
	os.system(cmd_samtools_index_1)
	print cmd_select_chr
	os.system(cmd_select_chr)
	print cmd_samtools_2
	os.system(cmd_samtools_2)
	print cmd_picard
	os.system(cmd_picard)
	print cmd_samtools_index_2
	os.system(cmd_samtools_index_2)
	print cmd_add_readgroup
	os.system(cmd_add_readgroup)
	print cmd_buildbamindex
	os.system(cmd_buildbamindex)
	print cmd_RealignerTargetCreator
	os.system(cmd_RealignerTargetCreator)
	print cmd_IndelRealigner
	os.system(cmd_IndelRealigner)
	print cmd_BaseRecalibrator
	os.system(cmd_BaseRecalibrator)
	print cmd_PrintReads
	os.system(cmd_PrintReads)
	

def netMHCpan(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir):
	str_proc=r'''
set -x
input_fasta=%s
hla_str=%s
netmhc_out=%s
out_dir=%s
split_num=%s
netMHCpan=%s
tmp=%s
if [ -d ${out_dir}/${tmp} ];then
	rm -rf ${out_dir}/${tmp}
	mkdir -p ${out_dir}/${tmp}
else
	mkdir -p ${out_dir}/${tmp}
fi
if [ -f ${netmhc_out} ];then
	rm ${netmhc_out}
fi
split -l ${split_num} ${input_fasta} ${out_dir}/${tmp}/
filelist=`ls ${out_dir}/${tmp}/`
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
		$netMHCpan -a $s -f ${out_dir}/${tmp}/${file_l} -l 8,9,10,11 > ${out_dir}/${tmp}/${s}_${file_l}_tmp_netmhc.txt
	} &
	done
	wait
}
done
for file_l in ${arr1[@]}
do
{
	rm ${out_dir}/${tmp}/${file_l}
}
done
filelist1=`ls ${out_dir}/${tmp}/`
for file_r in $filelist1
do
{
	cat ${out_dir}/${tmp}/${file_r} >> ${netmhc_out}
	rm ${out_dir}/${tmp}/${file_r}	
}
done
rm -rf 	${out_dir}/${tmp}
set +x
'''%(fasta_file,hla_str,netmhc_out_file,out_dir,split_num,netMHCpan_path,tmp_dir)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_somatic_caling_drift(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold):
	str_proc = r'''
set -e
somat_f=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
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
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_normal.bam  > ${PREFIX}_normal.fifo &
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_tumor.bam  > ${PREFIX}_tumor.fifo &
java -jar $varscan somatic ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX} #--output-vcf 1
java -jar $varscan processSomatic ${PREFIX}.snp
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
sed '1d' ${somat_f}/${PREFIX}.snp.Somatic | awk -F '\t' '{print $1,$2,$2,$3"/"$4}' > ${somat_f}/${PREFIX}_snv_vep_input.vcf
$vep -i ${somat_f}/${PREFIX}_snv_vep_input.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol --offline -o ${somat_f}/${PREFIX}_snv_all_vep_ann.txt
$vep -i ${somat_f}/${PREFIX}_snv_vep_input.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o ${somat_f}/${PREFIX}_snv_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/snv2fasta.py -i ${somat_f}/${PREFIX}_snv_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}
'''%(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_neo(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	#netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_snv')
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

def indel_calling_drift(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path):
	str_proc2=r'''
set -e
strelka_fold=%s
strelka_path=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
cpu=%s
vep=%s
if [ -d ${strelka_fold} ];then
	rm -rf ${strelka_fold}
fi
python ${strelka_path}/configureStrelkaSomaticWorkflow.py --tumorBam=${alignment_fold}/${PREFIX}_tumor_recal.bam --normalBam=${alignment_fold}/${PREFIX}_normal_recal.bam --referenceFasta=${REFERENCE} --config=${strelka_path}/configureStrelkaSomaticWorkflow.py.ini --runDir=${strelka_fold} --exome
python ${strelka_fold}/runWorkflow.py -m local -j $cpu -q ${PREFIX}_strelka -g 32 --quiet
$vep -i ${strelka_fold}/results/variants/somatic.indels.vcf.gz --cache --dir ${vep_cache} --dir_cache ${vep_cache} --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_strelka
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_strelka
'''%(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path)
	print str_proc2
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')
	
def indel_neo(somatic_mutation_fold,PREFIX,vep_cache,netmhc_out_fold,vep_path,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	str_proc1='''
somatic_mutation=%s
PREFIX=%s
vep_cache=%s
netmhc_out=%s
vep=%s
python ${iTuNES_BIN_PATH}/varscan_indel_preprocess.py -i ${somatic_mutation}/${PREFIX}.indel -o ${somatic_mutation} -s ${PREFIX}
$vep -i ${somatic_mutation}/${PREFIX}_varscan_indel.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(somatic_mutation_fold,PREFIX,vep_cache,netmhc_out_fold,vep_path)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')	
	str_proc3=r'''
PREFIX=%s
netmhc_out=%s
cat ${netmhc_out}/${PREFIX}_strelka_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_strelka_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_del.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
cat ${netmhc_out}/${PREFIX}_varscan_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
'''%(PREFIX,netmhc_out_fold)
	subprocess.call(str_proc3, shell=True, executable='/bin/bash')
	netMHCpan(indel_fasta_file,hla_str,indel_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_indel')
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
'''%(PREFIX,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	#print str_proc4
	subprocess.call(str_proc4, shell=True, executable='/bin/bash')
	
	
def varscan_copynumber_calling(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path):
	str_proc=r'''
set -x
copynumber_profile=%s
PREFIX=%s
alignments=%s
REFERENCE=%s
samtools=%s
varscan=%s
if [ ! -d ${copynumber_profile} ];then
	mkdir ${copynumber_profile}	
fi
rm -rf ${copynumber_profile}/*
cd ${copynumber_profile}
mkfifo ${PREFIX}_normal.fifo
mkfifo ${PREFIX}_tumor.fifo
$samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ${alignments}/${PREFIX}_normal_recal.bam > ${PREFIX}_normal.fifo &
$samtools mpileup -f ${REFERENCE} -q 5 -L 10000 -d 10000 ${alignments}/${PREFIX}_tumor_recal.bam > ${PREFIX}_tumor.fifo &
java -jar $varscan copynumber ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo ${PREFIX}
rm ${PREFIX}_normal.fifo ${PREFIX}_tumor.fifo
cd ..
'''%(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path)
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def pyclone_annotation(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_fold,coverage,pyclone_path,cancer_type):
	str_proc=r'''
somatic_mutation=%s
copynumber_profile=%s
PREFIX=%s
pyclone=%s
netctl=%s
COVERAGE=%d
Pyclone=%s
cancer_type=%s
#Rscript ${iTuNES_BIN_PATH}/sequenza_test.R ${somatic_mutation}/${PREFIX}.snp ${copynumber_profile}/${PREFIX}.copynumber ${copynumber_profile}/ ${PREFIX}
python ${iTuNES_BIN_PATH}/pyclone_input.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann.txt -s ${somatic_mutation}/${PREFIX}.snp.Somatic -c ${copynumber_profile}/${PREFIX}_seg_copynumber.txt -o ${pyclone} -S ${PREFIX} -C ${COVERAGE}
TUMOR_CONTENT=`cat ${copynumber_profile}/${PREFIX}_cellularity.txt`
$Pyclone run_analysis_pipeline --in_files ${pyclone}/${PREFIX}_pyclone_input.tsv --tumour_contents $TUMOR_CONTENT --prior major_copy_number --working_dir ${pyclone}
python ${iTuNES_BIN_PATH}/neo_pyclone_annotation.py -n ${netctl}/${PREFIX}_snv_netctl_concact.txt -i ${somatic_mutation}/${PREFIX}_snv_vep_ann.txt -s ${pyclone}/tables/loci.tsv -o ${netctl} -S ${PREFIX} -t ${cancer_type}
'''%(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_fold,coverage,pyclone_path,cancer_type)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')


def kallisto_expression(raw_fastq_path_first,raw_fastq_path_second,kallisto_path,kallisto_out_fold,prefix,kallisto_index_path):
	cmd_kallisto = kallisto_path + " quant -i " + kallisto_index_path + " -o " + kallisto_out_fold + " " + raw_fastq_path_first + " " + raw_fastq_path_second
	print cmd_kallisto
	os.system(cmd_kallisto)



def immunogenicity_score_calculate(final_neo_file,immunogenicity_score_ranking_file):
	data_neo = pd.read_table(final_neo_file,header=0,sep='\t')
	neoantigen_infor = []
	for i in range(len(data_neo.Gene)):
		neo_infor = data_neo["Gene"][i]+'_'+data_neo["AA_change"][i]+'_'+data_neo["MT_pep"][i]+'_'+data_neo["WT_pep"][i]
		neoantigen_infor.append(neo_infor)
	data_neo["neoantigen_infor"] = neoantigen_infor
	f_affinity_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(x-2))))/2
	f_affinity_rank_mt=lambda x:1/(1+math.pow(math.e,5*(x-2)))
	aff_mt_rank_score=data_neo.MT_Binding_level.apply(f_affinity_rank_mt)
	aff_wt_rank_score=data_neo.WT_Binding_level.apply(f_affinity_rank_wt)
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	tpm_score=data_neo.tpm.apply(f_TPM)
	f_normal_TPM=lambda x:1-math.tanh(x/k)
	tpm_normal_score=data_neo.iloc[:,18].apply(f_normal_TPM)
	allele_frequency_score=data_neo.variant_allele_frequency
	netchop_score=data_neo.combined_prediction_score
	cellular_prevalence_score=data_neo.cellular_prevalence
	immunogenicity_score=[]
	for i in range(len(aff_mt_rank_score)):
		IS=aff_mt_rank_score[i]*aff_wt_rank_score[i]*tpm_score[i]*tpm_normal_score[i]*allele_frequency_score[i]*netchop_score[i]*cellular_prevalence_score[i]
		immunogenicity_score.append(IS)
	data_neo["immunogenicity_score"]=immunogenicity_score
	data_sort=data_neo.sort_values(["immunogenicity_score"],ascending=False)
	data_sort.to_csv(immunogenicity_score_ranking_file,header=1,index=0,sep='\t')

def rankagg_score_calculate(final_neo_file,neo_for_rankagg_file,neo_rank_agg_result):
	data_neo = pd.read_table(final_neo_file,header=0,sep='\t')
	neoantigen_infor = []
	for i in range(len(data_neo.Gene)):
		neo_infor = data_neo["Gene"][i]+'_'+data_neo["AA_change"][i]+'_'+data_neo["MT_pep"][i]+'_'+data_neo["WT_pep"][i]
		neoantigen_infor.append(neo_infor)
	data_neo["neoantigen_infor"] = neoantigen_infor
	f_affinity_rank_wt=lambda x:1-(1/(1+math.pow(math.e,5*(x-2))))/2
	f_affinity_rank_mt=lambda x:1/(1+math.pow(math.e,5*(x-2)))
	aff_mt_rank_score=data_neo.MT_Binding_level.apply(f_affinity_rank_mt)
	aff_wt_rank_score=data_neo.WT_Binding_level.apply(f_affinity_rank_wt)
	k=1
	f_TPM=lambda x:math.tanh(x/k)
	tpm_score=data_neo.tpm.apply(f_TPM)
	f_normal_TPM=lambda x:1-math.tanh(x/k)
	tpm_normal_score=data_neo.iloc[:,18].apply(f_normal_TPM)
	allele_frequency_score=data_neo.variant_allele_frequency
	netchop_score=data_neo.combined_prediction_score
	cellular_prevalence_score=data_neo.cellular_prevalence
	data_rank_agg=pd.DataFrame()
	data_rank_agg["neoantigen_infor"]=neoantigen_infor
	data_rank_agg["aff_mt_rank_score"]=aff_mt_rank_score
	data_rank_agg["aff_wt_rank_score"]=aff_wt_rank_score
	data_rank_agg["tpm_score"]=tpm_score
	data_rank_agg["tpm_normal_score"]=tpm_normal_score
	data_rank_agg["allele_frequency_score"]=allele_frequency_score
	data_rank_agg["netchop_score"]=netchop_score
	data_rank_agg["cellular_prevalence_score"]=cellular_prevalence_score
	data_rank_agg.to_csv(neo_for_rankagg_file,header=1,sep='\t',index=0)
	r=R()
	str_proc='''
neo_for_rankagg=\"%s\"
neo_rank_agg=\"%s\"
library(tidyr)
library("RobustRankAggreg")
data=read.table(neo_for_rankagg,sep='\t',head=TRUE)
list_mt_binding=as.vector(data[order(data$aff_mt_rank_score,decreasing=T),]['neoantigen_infor'][[1]])
list_wt_binding=as.vector(data[order(data$aff_wt_rank_score,decreasing=T),]['neoantigen_infor'][[1]])
list_TPM=as.vector(data[order(data$tpm_score,decreasing=T),]['neoantigen_infor'][[1]])
list_combined_prediction_score=as.vector(data[order(data$netchop_score,decreasing=T),]['neoantigen_infor'][[1]])
list_cellular_prevalence=as.vector(data[order(data$cellular_prevalence_score,decreasing=T),]['neoantigen_infor'][[1]])
list_variant_allele_frequency=as.vector(data[order(data$allele_frequency_score,decreasing=T),]['neoantigen_infor'][[1]])
list_normal_expression=as.vector(data[order(data$tpm_normal_score,decreasing=T),]['neoantigen_infor'][[1]])
glist=list(list_mt_binding,list_wt_binding,list_TPM,list_combined_prediction_score,list_cellular_prevalence,list_variant_allele_frequency,list_normal_expression)
r=rankMatrix(glist)
agg_result=aggregateRanks(rmat = r, method = "RRA")
colnames(agg_result)=c("Ranked_neoantigen_name","Ranking_score")
write.table(agg_result,file=neo_rank_agg,row.names=F,quote=F,sep='\t')
'''%(neo_for_rankagg_file,neo_rank_agg_result)
	r(str_proc)

