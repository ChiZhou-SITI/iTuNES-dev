import os
from VCFprocessor import *
import shutil
def Vcf(opts):
	vcf_file=opts.vcf
	hla_type=opts.hlatype
	out_directory=opts.dir
	prefix=opts.PREFIX
	cpu_thread=opts.CPU
	binding_affinity_cutoff=opts.B_CUTOFF
	binding_aff_fc_cutoff=opts.A_CUTOFF
	exp=opts.EXPRESSION_FILE
	fpkm_cutoff=opts.FPKM_CUTOFF
	neoantigens_vcf(vcf_file, hla_type ,out_directory,prefix,cpu_thread,binding_aff_fc_cutoff,binding_affinity_cutoff,exp,fpkm_cutoff)	
print 'Done!'

#	if os.path.exists(str_path_mutsigcv_aux):
#		shutil.rmtree(str_path_mutsigcv_aux)


	
