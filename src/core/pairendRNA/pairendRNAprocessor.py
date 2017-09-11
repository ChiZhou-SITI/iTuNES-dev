import os,sys,time
import multiprocessing
import shutil 
import subprocess


def hlatyping(raw_fastq_path_first,raw_fastq_path_second,opitype_fold,opitype_out_fold,opitype_ext,prefix):
	os.system("rm -rf %s/*"%opitype_out_fold)
	cmd_hla = 'python ' + opitype_fold + ' -i ' + raw_fastq_path_first + ' ' + raw_fastq_path_second + ' --rna -o ' + opitype_out_fold
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

def mapping_qc_gatk_preprocess(fastq_1_path,fastq_2_path,CPU,STAR_index,alignment_out_fold,prefix,REFERENCE,STAR_path,stringtie_path,java_picard_path,GATK_path,dbsnp138,OneKG,mills,GTF_path,expression_fold):
	cmd_STAR=STAR_path + " --runThreadN " + str(CPU) + " --twopassMode Basic --readFilesIn " + fastq_1_path + " " +  fastq_2_path + " --genomeDir " + STAR_index + " --outFileNamePrefix " + alignment_out_fold+'/'+prefix + " --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64424509440 "
	cmd_stringtie=stringtie_path + alignment_out_fold+'/'+prefix+"Aligned.sortedByCoord.out.bam" + " -o " + expression_fold+'/'+prefix+".gtf" + " -p 8 -G " + GTF_path + " -A " + expression_fold + '/' + prefix + "_gene_abund.tab"
	cmd_add_readgroup="java -Xmx4G -jar " + java_picard_path + ' AddOrReplaceReadGroups I=' + alignment_out_fold+'/'+ prefix + 'Aligned.sortedByCoord.out.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_add.bam' + ' SO=coordinate VALIDATION_STRINGENCY=SILENT RGID=id RGLB=solexa-123 RGPL=illumina RGPU=AXL2342  RGSM=WGC015802 RGCN=bi RGDT=2014-01-20'
	cmd_picard="java -Xmx4G -jar " + java_picard_path + ' MarkDuplicates INPUT=' + alignment_out_fold+'/'+ prefix + '_add.bam' + ' OUTPUT=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' METRICS_FILE=' + alignment_out_fold+'/'+prefix + '_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'
	cmd_buildbamindex="java -Xmx4G -jar " + java_picard_path + ' BuildBamIndex I=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' O=' + alignment_out_fold+'/'+ prefix + '_add_markdup.bam.bai' + ' VALIDATION_STRINGENCY=SILENT'
	cmd_SplitNCigarReads="java -Xmx4G -jar " + GATK_path + ' -T SplitNCigarReads -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_add_markdup.bam' + ' -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o ' + alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam'
	cmd_RealignerTargetCreator="java -Xmx4G -jar " + GATK_path + ' -T RealignerTargetCreator -nt 8 -dt NONE -R ' + REFERENCE + ' -I '+ alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '.intervals'
	cmd_IndelRealigner="java -Xmx4G -jar " + GATK_path + ' -T IndelRealigner -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_SplitNCigarReads.bam' + ' -known ' + OneKG + ' -known ' + mills + ' -targetIntervals ' + alignment_out_fold+'/'+ prefix + '.intervals' + ' -filterRNC -o ' + alignment_out_fold+'/'+ prefix + '_realign.bam'
	cmd_BaseRecalibrator="java -Xmx4G -jar " + GATK_path + ' -T BaseRecalibrator -nct 8 -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix +'_realign.bam' + ' -knownSites ' + OneKG + ' -knownSites ' + mills + ' -knownSites ' + dbsnp138 + ' -o ' + alignment_out_fold + '/' + prefix + '.table'
	cmd_PrintReads="java -Xmx4G -jar " + GATK_path + ' -T PrintReads -nct 8 -dt NONE -R ' + REFERENCE + ' -I ' + alignment_out_fold+'/'+ prefix + '_realign.bam' + ' -BQSR ' + alignment_out_fold + '/' + prefix + '.table' + ' -o ' + alignment_out_fold + '/' + prefix + '_recal.bam'
	print cmd_STAR
	os.system(cmd_STAR)
	print cmd_stringtie
	os.system(cmd_stringtie)
	print cmd_add_readgroup
	os.system(cmd_add_readgroup)
	print cmd_picard
	os.system(cmd_picard)
	print cmd_buildbamindex
	os.system(cmd_buildbamindex)
	print cmd_SplitNCigarReads
	os.system(cmd_SplitNCigarReads)
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


def varscan_snv_calling(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold):
	str_proc = r'''
set -e
somatic_mutation=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
if [ ! -d ${somatic_mutation} ];then
	mkdir ${somatic_mutation}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
rm -rf ${somatic_mutation}/*
cd ${somatic_mutation}
mkfifo ${PREFIX}_tumor.fifo
$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_recal.bam > ${PREFIX}_tumor.fifo &
java -jar $varscan mpileup2snp ${PREFIX}_tumor.fifo --output-vcf 1 > ${somatic_mutation}/${PREFIX}_varscan_snv.vcf #--output-vcf 1
rm ${PREFIX}_tumor.fifo
cd ..
$vep -i ${somatic_mutation}/${PREFIX}_varscan_snv.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is missense_variant" -o ${somatic_mutation}/${PREFIX}_varscan_snv_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/snv2fasta.py -i ${somatic_mutation}/${PREFIX}_varscan_snv_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(somatic_mutation_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')
def varscan_indel_calling(varscan_indel_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold):
	str_proc = r'''
set -e
varscan_indel=%s
alignment_fold=%s
PREFIX=%s
REFERENCE=%s
vep_cache=%s
netmhc_out=%s
samtools=%s
varscan=%s
vep=%s
if [ ! -d ${varscan_indel} ];then
	mkdir ${varscan_indel}	
fi
if [ ! -d ${netmhc_out} ];then
	mkdir ${netmhc_out}	
fi
#rm -rf ${varscan_indel}/*
#cd ${varscan_indel}
#mkfifo ${PREFIX}_tumor.fifo
#$samtools mpileup -f ${REFERENCE} -q 5 -Q 20 -L 10000 -d 10000 ${alignment_fold}/${PREFIX}_recal.bam  > ${PREFIX}_tumor.fifo &
#java -jar $varscan mpileup2indel ${PREFIX}_tumor.fifo --output-vcf 1 > ${varscan_indel}/${PREFIX}_varscan_indel.vcf #--output-vcf 1
#rm ${PREFIX}_tumor.fifo
#cd ..
$vep -i ${varscan_indel}/${PREFIX}_varscan_indel.vcf --cache --dir $vep_cache --dir_cache $vep_cache --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_varscan
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${varscan_indel}/${PREFIX}_varscan_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_varscan
'''%(varscan_indel_fold,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,samtools_path,varscan_path,vep_path)
	print str_proc
	subprocess.call(str_proc, shell=True, executable='/bin/bash')

def strelka_indel_calling(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path):
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
python ${strelka_path}/configureStrelkaGermlineWorkflow.py --bam=${alignment_fold}/${PREFIX}_recal.bam --referenceFasta=${REFERENCE} --config=${strelka_path}/configureStrelkaGermlineWorkflow.py.ini --runDir=${strelka_fold} --rna
python ${strelka_fold}/runWorkflow.py -m local -j $cpu -q ${PREFIX}_strelka -g 32 --quiet
$vep -i ${strelka_fold}/results/variants/germline.indels.vcf.gz --cache --dir ${vep_cache} --dir_cache ${vep_cache} --force_overwrite  --symbol -o STDOUT --offline | filter_vep --ontology --filter "Consequence is coding_sequence_variant" -o ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt --force_overwrite
python ${iTuNES_BIN_PATH}/varscandel2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt -o ${netmhc_out} -s ${PREFIX}_strelka
python ${iTuNES_BIN_PATH}/varscanins2fasta.py -i ${strelka_fold}/${PREFIX}_strelka_indel_vep_ann.txt  -o ${netmhc_out} -s ${PREFIX}_strelka
'''%(strelka_out_fold,strelka_path,alignment_out_fold,PREFIX,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path)
	print str_proc2
	subprocess.call(str_proc2, shell=True, executable='/bin/bash')
	
def varscan_neo(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	netMHCpan(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_snv')
	str_proc1=r'''
PREFIX=%s
netmhc_out=%s
Exp_file=%s
Binding_Aff_Fc_Cutoff=%d
Binding_Aff_Cutoff=%d
Fpkm_Cutoff=%d
hla_str=%s
netctl_fold=%s
python ${iTuNES_BIN_PATH}/sm_netMHC_result_parse.py -i ${netmhc_out}/${PREFIX}_snv_netmhc.txt -g ${netmhc_out}/${PREFIX}_varscan_snv.fasta -o ${netmhc_out} -s ${PREFIX}_snv -e ${Exp_file} -a ${Binding_Aff_Fc_Cutoff} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/netCTLPAN.py -i ${netmhc_out}/${PREFIX}_snv_final_neo_candidate.txt -o ${netctl_fold} -s ${PREFIX}_snv
'''%(prefix,netmhc_out_fold,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,hla_str,netctl_fold)
	print str_proc1
	subprocess.call(str_proc1, shell=True, executable='/bin/bash')

def indel_neo(PREFIX,netmhc_out_fold,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_fold,netMHCpan_path):
	str_proc3=r'''
PREFIX=%s
netmhc_out=%s
#cat ${netmhc_out}/${PREFIX}_strelka_del.fasta > ${netmhc_out}/${PREFIX}_indel.fasta
#cat ${netmhc_out}/${PREFIX}_strelka_ins.fasta >> ${netmhc_out}/${PREFIX}_indel.fasta
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

def gene_fusion(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold,genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,split_num,netMHCpan_path,Binding_Aff_Cutoff,Fpkm_Cutoff):
	str_proc='''
T_fastq_1=%s
T_fastq_2=%s
prefix=%s
gene_fusion_fold=%s
eric_db=%s
eric_path=%s
cpu=%d
netmhc_out=%s
$eric_path -db $eric_db --refid homo_sapiens -name ${prefix}_genefusion -o result $T_fastq_1 $T_fastq_2 -p $cpu
cd ..
python ${iTuNES_BIN_PATH}/genefusion2fasta.py -i ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -c ${GeneFusion_Cutoff} -o ${netmhc_out} -s ${prefix}
'''%(fastq_1_path,fastq_2_path,prefix,gene_fusion_fold_path,eric_db_path,eric_path,cpu,netmhc_out_fold)
	netMHCpan(genefuison_fasta_file,hla_str,genefuison_netmhc_out_file,netmhc_out_fold,split_num,netMHCpan_path,'tmp_genefuison')
	str_proc1='''
netmhc_out=%s
prefix=%s
Binding_Aff_Cutoff=%s
Fpkm_Cutoff=%s
hla_str=%s
python ${iTuNES_BIN_PATH}/gf_netMHC_result_parse.py -i ${netmhc_out}/${prefix}_genefusion_netmhc.txt -t ${netmhc_out}/${prefix}_gene_fusion.fasta -g ${gene_fusion_fold}/result/${prefix}_genefusion.results.total.tsv -o ${netmhc_out} -s ${prefix} -b ${Binding_Aff_Cutoff} -f ${Fpkm_Cutoff} -l ${hla_str}
python ${iTuNES_BIN_PATH}/gf_netctl.py -i netmhc/${prefix}_genefusion_final_neo_candidate.txt -o netctl -s ${prefix}_gf
'''%(netmhc_out_fold,prefix,Binding_Aff_Cutoff,Fpkm_Cutoff,hla_str)

