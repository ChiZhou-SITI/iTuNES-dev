import argparse

def ParsePAIRRNA(p_perna):
	p_perna_inp = p_perna.add_argument_group('Input options')
	p_perna_inp.add_argument('-p', '--input_tumor_fastq1', dest='tumor_fq1', required=True, help='Input first fastq file from tumor sample(required)')
	p_perna_inp.add_argument('-q', '--input_tumor_fastq2', dest='tumor_fq2', required=True, help='Input secend fastq file from tumor sample(required)')
	p_perna_inp.add_argument('-o', '--output_dir', dest='dir', required=True,help='Output all result in this fold.')
	p_perna_otp = p_perna.add_argument_group('Optional arguments')
	p_perna_otp.add_argument('-g', '--genome', dest='GENOME', help='specify the species,defalut:human')
	p_perna_otp.add_argument('-G', '--gtf', dest='GTF', help='specify the annotation gtf file,default=/home/zhouchi/database/Annotation/GTF/human.gtf')
	p_perna_otp.add_argument('-P', '--prefix', dest='PREFIX', help='specify the prefix of output file,default: input fastq file name')
	p_perna_otp.add_argument('-r', '--reference', dest='REFERENCE', help='specify the reference fasta file of specified species,default=/home/zhouchi/database/Annotation/Fasta/human.fasta')
	p_perna_otp.add_argument('-I', '--star_index', dest='STAR_INDEX', help='specify the path of STAR index,default: /home/zhouchi/database/Annotation/Index/STAR/human')
	p_perna_otp.add_argument('-c', '--cpu', dest='CPU', help='specify the thread,default=4',default=4)
	p_perna_otp.add_argument('-C', '--coverage', dest='COVERAGE', help='specify the coverage cutoff of sequence,default=200',default=200)
	p_perna_otp.add_argument('-a', '--bindingaffinity_foldchange_cutoff', dest='A_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=1)
	p_perna_otp.add_argument('-b', '--bindingaffinity_cutoff', dest='B_CUTOFF', help='specify the binding affinity cutoff of candidate neoantigen,default=500nM',default=500)
	p_perna_otp.add_argument('-E', '--expression_file', dest='EXPRESSION_FILE', help='the expression profile with FPKM value,default=500nM')
	p_perna_otp.add_argument('-f', '--FPKM_cutoff', dest='FPKM_CUTOFF', help='set FPKM value cutoff to filter lowly expressed neoantigens,set this option together with -E,default=2',default=2)

	
	
