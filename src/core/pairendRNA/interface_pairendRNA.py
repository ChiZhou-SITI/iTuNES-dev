import os
from pairendRNAprocessor import *
import shutil
def PERNA(opts):
	tumor_fastq_1=opts.tumor_fq1
	tumor_fastq_2=opts.tumor_fq2
	out_directory=opts.dir
	genome=opts.GENOME
	prefix=opts.PREFIX
	gtf_file=opts.GTF
	reference=opts.REFERENCE
	star_index=opts.STAR_INDEX
	cpu_thread=opts.CPU
	coverage=opts.COVERAGE
	binding_affinity_cutoff=opts.B_CUTOFF
	binding_aff_fc_cutoff=opts.A_CUTOFF
	exp=opts.EXPRESSION_FILE
	fpkm_cutoff=opts.FPKM_CUTOFF
	print binding_affinity_cutoff,exp,binding_aff_fc_cutoff
	neoantigens_PairRna(tumor_fastq_1,tumor_fastq_2,out_directory,genome,prefix,gtf_file,star_index,reference,cpu_thread,coverage,binding_aff_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff)	
print 'Done!'

#	if os.path.exists(str_path_mutsigcv_aux):
#		shutil.rmtree(str_path_mutsigcv_aux)


	
