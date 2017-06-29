import subprocess
import os,sys,time

def neoantigens_SingleNoMatchDna(tumor_fastq,out_directory,genome,prefix,bwa_index,reference,cpu_thread,coverage,binding_affinity_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff):
	str_proc = '''
T_fq=%s
out_dir=%s
genome_type=%s
pre=%s
index_bwa=%s
reference_fasta=%s
cpu_num=%s
coverage_cutoff=%s
ba_fc_cutoff=%s
ba_cutoff=%s
exp_file=%s
fpkm_cutoff=%s
bash ${iTuNES_BIN_PATH}/neoantigen_pipeline_DNAnoMatched_se.sh \
-p ${T_fq} \
-o ${out_dir} \
-g ${genome_type} \
-P ${pre} \
-I ${index_bwa} \
-r ${reference_fasta} \
-c ${cpu_num} \
-C ${coverage_cutoff} \
-a ${ba_fc_cutoff} \
-b ${ba_cutoff} \
-E ${exp_file} \
-f ${fpkm_cutoff}
	'''%(tumor_fastq,out_directory,genome,prefix,bwa_index,reference,cpu_thread,coverage,binding_affinity_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff)
	subprocess.call(str_proc,shell=True,executable='/bin/bash')














