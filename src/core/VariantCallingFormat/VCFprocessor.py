import subprocess
import os,sys,time

def neoantigens_vcf(input_vcf_file,hla_type,out_directory,prefix,cpu_thread,binding_affinity_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff):
	str_proc = '''
vcf_in=%s
hlatype=%s
out_dir=%s
pre=%s
cpu_num=%s
ba_fc_cutoff=%s
ba_cutoff=%s
exp_file=%s
fpkm_cutoff=%s
bash ${iTuNES_BIN_PATH}/neoantigen_pipeline_vcf.sh \
-i ${vcf_in} \
-t ${hlatype} \
-o ${out_dir} \
-P ${pre} \
-c ${cpu_num} \
-a ${ba_fc_cutoff} \
-b ${ba_cutoff} \
-E ${exp_file} \
-f ${fpkm_cutoff}
	'''%(input_vcf_file,hla_type,out_directory,prefix,cpu_thread,binding_affinity_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff)
	subprocess.call(str_proc,shell=True,executable='/bin/bash')














