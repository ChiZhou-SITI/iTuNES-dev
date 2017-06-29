import argparse

def ParseVCF(p_vcf):
	p_vcf_inp = p_vcf.add_argument_group('Input options')
	p_vcf_inp.add_argument('-i', '--input_vcf', dest='vcf', required=True, help='Input mutation data in vcf format')
	p_vcf_inp.add_argument('-t', '--hla_type', dest='hlatype', required=True, help='Input hla type of the sample')
	p_vcf_inp.add_argument('-o', '--output_dir', dest='dir', required=True,help='Output all result in this fold.')
	p_vcf_otp = p_vcf.add_argument_group('Optional arguments')
	p_vcf_otp.add_argument('-P', '--prefix', dest='PREFIX', help='specify the prefix of output file,default: input fastq file name')
	p_vcf_otp.add_argument('-c', '--cpu', dest='CPU', help='specify the thread,default=4',default=4)
	p_vcf_otp.add_argument('-a', '--bindingaffinity_foldchange_cutoff', dest='A_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=1)
	p_vcf_otp.add_argument('-b', '--bindingaffinity_cutoff', dest='B_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=500)
	p_vcf_otp.add_argument('-E', '--expression_file', dest='EXPRESSION_FILE', help='the expression profile with FPKM value,default=500nM')
	p_vcf_otp.add_argument('-f', '--FPKM_cutoff', dest='FPKM_CUTOFF', help='set FPKM value cutoff to filter lowly expressed neoantigens,set this option together with -E,default=2',default=2)