import os
from singleendNMDNAprocessor import *
import shutil
def SENMD(opts):
	tumor_fastq=opts.tumor_fq
	out_directory=opts.dir
	genome=opts.GENOME
	prefix=opts.PREFIX
	reference=opts.REFERENCE
	bwa_index=opts.BWA_INDEX
	cpu_thread=opts.CPU
	coverage=opts.COVERAGE
	binding_affinity_cutoff=opts.B_CUTOFF
	binding_aff_fc_cutoff=opts.A_CUTOFF
	exp=opts.EXPRESSION_FILE
	fpkm_cutoff=opts.FPKM_CUTOFF
	neoantigens_SingleNoMatchDna(tumor_fastq,out_directory,genome,prefix,bwa_index,reference,cpu_thread,coverage,binding_aff_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff)	

#	if os.path.exists(str_path_mutsigcv_aux):
#		shutil.rmtree(str_path_mutsigcv_aux)


	
