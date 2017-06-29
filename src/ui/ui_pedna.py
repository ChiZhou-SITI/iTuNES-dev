import argparse

def ParsePAIRMATCHDNA(p_pedna):
	p_pedna_inp = p_pedna.add_argument_group('Input options')
	p_pedna_inp.add_argument('-p', '--input_normal_fastq1', dest='normal_fq1', required=True, help='Input first fastq file from normal sample(required)')
	p_pedna_inp.add_argument('-q', '--input_normal_fastq2', dest='normal_fq2', required=True, help='Input secend fastq file from normal sample(required)')
	p_pedna_inp.add_argument('-s', '--input_tumor_fastq1', dest='tumor_fq1', required=True, help='Input first fastq file from tumor sample(required)')
	p_pedna_inp.add_argument('-t', '--input_tumor_fastq2', dest='tumor_fq2', required=True, help='Input secend fastq file from normal sample(required)')
	p_pedna_inp.add_argument('-o', '--output_dir', dest='dir', required=True,help='Output all result in this fold.')
	p_pedna_otp = p_pedna.add_argument_group('Optional arguments')
	p_pedna_otp.add_argument('-g', '--genome', dest='GENOME', help='specify the species,defalut:human')
	p_pedna_otp.add_argument('-P', '--prefix', dest='PREFIX', help='specify the prefix of output file,default: input fastq file name')
	p_pedna_otp.add_argument('-r', '--reference', dest='REFERENCE', help='specify the reference fasta file of specified species,default=/home/zhouchi/database/Annotation/Fasta/human.fasta')
	p_pedna_otp.add_argument('-I', '--bwa_index', dest='BWA_INDEX', help='specify the path of bwa index,default: /home/zhouchi/database/Annotation/Index/bwa/human')
	p_pedna_otp.add_argument('-c', '--cpu', dest='CPU', help='specify the thread,default=4',default=4)
	p_pedna_otp.add_argument('-C', '--coverage', dest='COVERAGE', help='specify the coverage cutoff of sequence,default=200',default=5)
	p_pedna_otp.add_argument('-a', '--bindingaffinity_foldchange_cutoff', dest='A_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=1)
	p_pedna_otp.add_argument('-b', '--bindingaffinity_cutoff', dest='B_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=500)
	p_pedna_otp.add_argument('-E', '--expression_file', dest='EXPRESSION_FILE', help='the expression profile with FPKM value,default=500nM')
	p_pedna_otp.add_argument('-f', '--FPKM_cutoff', dest='FPKM_CUTOFF', help='set FPKM value cutoff to filter lowly expressed neoantigens,set this option together with -E,default=2',default=2)

	
	
