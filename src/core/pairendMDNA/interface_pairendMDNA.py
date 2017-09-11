import os
from pairendMDNAprocessor import *
import shutil
import yaml
def PEMD(opts):
	config_file=opts.Config_file
	f=open(config_file)
	config_list=yaml.load(f)
	#######read and parse parameter
	print "read and parse parameter."
	output_fold=config_list["output_fold"]
	itunes_bin_path=config_list["itunes_bin_path"]
	os.system("export iTuNES_BIN_PATH=%s"%itunes_bin_path)
	normal_fastq_path_first=config_list["normal_fastq_path_first"]
	normal_fastq_path_second=config_list["normal_fastq_path_second"]
	tumor_fastq_path_first=config_list["tumor_fastq_path_first"]
	tumor_fastq_path_second=config_list["tumor_fastq_path_second"]
	opitype_fold=config_list["opitype_fold"]
	opitype_out_fold=output_fold + '/' + 'hlatyping'
	opitype_ext='${iTuNES_BIN_PATH}/optitype_ext.py'
	prefix=config_list["sample_name"]
	CPU=config_list["thread_number"]
	REFERENCE=config_list["reference_path"]
	BWA_INDEX=config_list["bwa_index_path"]
	GENOME=config_list["genome"]
	alignment_out_fold=output_fold + '/' + 'alignments'
	bwa_path=config_list["bwa_path"]
	samtools_path=config_list["samtools_path"]
	java_picard_path=config_list["java_picard_path"]
	GATK_path=config_list["GATK_path"]
	dbsnp138_path=config_list["dbsnp138"]
	hapmap_path=config_list["hapmap"]
	omni_path=config_list["omni"]
	OneKG_path=config_list["OneKG"]
	mills_path=config_list["mills"]
	cosmic_path=config_list["cosmic"]
	somatic_mutation_fold=output_fold + '/' + 'somatic_mutation'
	vep_cache=config_list["vep_cache_path"]
	varscan_path=config_list["varscan_path"]
	vep_path=config_list["vep_path"]
	netmhc_out_fold=output_fold + '/' + 'netmhc'
	strelka_out_fold=output_fold + '/' + 'strelka'
	indel_fasta_file=netmhc_out_fold+'/'+prefix+'_indel.fasta'
	hla_str=config_list["hla_str"]
	indel_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_indel_netmhc.txt'
	split_num=1000
	binding_fc_aff_cutoff=int(config_list["binding_fc_aff_cutoff"])
	binding_aff_cutoff=int(config_list["binding_aff_cutoff"])
	fpkm_cutoff=int(config_list["fpkm_cutoff"])
	netctl_out_fold=output_fold + '/' + 'netctl'
	netMHCpan_path=config_list["netMHCpan_path"]
	varscan_copynumber_fold=output_fold + '/' + 'copynumber_profile'
	snv_fasta_file=netmhc_out_fold+'/'+prefix+'_snv.fasta'
	snv_netmhc_out_file=netmhc_out_fold+'/'+prefix+'_snv_netmhc.txt'
	pyclone_fold=output_fold + '/' + 'pyclone'
	coverage=int(config_list["coverage_cutoff"])
	pyclone_path=config_list["pyclone_path"]
	strelka_path=config_list["strelka_path"]
	cancer_type=config_list["cancer_type"]
	rna_fastq_1_path=config_list["tumor_rna_fastq_1"]
	rna_fastq_2_path=config_list["tumor_rna_fastq_2"]
	kallisto_path=config_list["kallisto_path"]
	kallisto_out_fold=output_fold + '/' + 'expression'
	kallisto_index_path=config_list["kallisto_index_path"]
	final_neo_file=netctl_out_fold + '/' + prefix + '_pyclone_neo.txt'
	candidate_neoantigens_fold=output_fold + '/' + 'candidate_neoantigens'
	immunogenicity_score_ranking_file=candidate_neoantigens_fold + '/' + prefix + '_immunogenicity_score_ranking.txt'
	neo_for_rankagg_file=netctl_out_fold + '/' + prefix + '_neo_for_rankagg.txt'
	neo_rank_agg_result=candidate_neoantigens_fold + '/'+ prefix + '_neo_rankagg_score_ranking.txt'
	trimmomatic_path=config_list["trimmomatic_path"]
	adapter_path=config_list["adapter_path"]
	tumor_fastq_prefix_first=os.path.splitext(os.path.basename(tumor_fastq_path_first))[0]
	tumor_fastq_prefix_second=os.path.splitext(os.path.basename(tumor_fastq_path_second))[0]
	tumor_fastq_dir=os.path.dirname(tumor_fastq_path_first)
	tumor_fastq_clean_first=fastq_dir + '/' + tumor_fastq_prefix_first + "_clean.fastq"
	tumor_fastq_clean_second=fastq_dir + '/' + tumor_fastq_prefix_second + "_clean.fastq"
	tumor_fastq_unpaired_first=fastq_dir + '/' + tumor_fastq_prefix_first + "_unpaired.fastq"
	tumor_fastq_unpaired_second=fastq_dir + '/' + tumor_fastq_prefix_second + "_unpaired.fastq"
	
	normal_fastq_prefix_first=os.path.splitext(os.path.basename(normal_fastq_path_first))[0]
	normal_fastq_prefix_second=os.path.splitext(os.path.basename(normal_fastq_path_second))[0]
	normal_fastq_dir=os.path.dirname(normal_fastq_path_first)
	normal_fastq_clean_first=fastq_dir + '/' + normal_fastq_prefix_first + "_clean.fastq"
	normal_fastq_clean_second=fastq_dir + '/' + normal_fastq_prefix_second + "_clean.fastq"
	normal_fastq_unpaired_first=fastq_dir + '/' + normal_fastq_prefix_first + "_unpaired.fastq"
	normal_fastq_unpaired_second=fastq_dir + '/' + normal_fastq_prefix_second + "_unpaired.fastq"


	print cancer_type
	if os.path.exists(rna_fastq_1_path) and os.path.exists(rna_fastq_2_path):
		exp_file=output_fold + '/' + "expression/abundance.tsv"
	else:
		exp_file=config_list["expression_file"]
	#####check input file,tool path and reference file#####
	if os.path.exists(normal_fastq_path_first) and os.path.exists(normal_fastq_path_second) and os.path.exists(tumor_fastq_path_first) and os.path.exists(tumor_fastq_path_second):
		print "check all fastq file done."
	else:
		print "please check your input fastq file!"
		os._exit(1)
	if os.path.exists(opitype_fold):
		print "check opitype path done."
	else:
		print "please check your opitype path!"
		os._exit(1)
	if os.path.exists(bwa_path):
		print "check bwa path done."
	else:
		print "please check your bwa path!"
		os._exit(1)	
	if os.path.exists(kallisto_path):
		print "check kallisto path done."
	else:
		print "please check your kallisto path!"
		os._exit(1)	
	if os.path.exists(samtools_path):
		print "check samtools path done."
	else:
		print "please check your samtools path!"
		os._exit(1)	
	if os.path.exists(java_picard_path):
		print "check picard path done."
	else:
		print "please check your picard path!"
		os._exit(1)	
	if os.path.exists(GATK_path):
		print "check GATK path done."
	else:
		print "please check your GATK path!"
		os._exit(1)	
	if os.path.exists(dbsnp138_path):
		print "check dbsnp138 file path done."
	else:
		print "please check your dbsnp138 file path!"
		os._exit(1)	
	if os.path.exists(OneKG_path):
		print "check OneKG file path done."
	else:
		print "please check your OneKG file path!"
		os._exit(1)	
	if os.path.exists(mills_path):
		print "check mills file path done."
	else:
		print "please check your mills file path!"
		os._exit(1)	
	if os.path.exists(vep_path):
		print "check vep path done."
	else:
		print "please check your vep path!"
		os._exit(1)	
	if os.path.exists(vep_cache):
		print "check vep cache path done."
	else:
		print "please check your vep cache path!"
		os._exit(1)		
	if os.path.exists(REFERENCE):
		print "check REFERENCE file path done."
	else:
		print "please check your REFERENCE file path!"
		os._exit(1)	
	if os.path.exists(BWA_INDEX+'/'+GENOME+".amb") and os.path.exists(BWA_INDEX+'/'+GENOME+".bwt") and os.path.exists(BWA_INDEX+'/'+GENOME+".pac") and os.path.exists(BWA_INDEX+'/'+GENOME+".sa"): 
		print "bwa index path:%s"%(BWA_INDEX)
	elif os.path.exists(REFERENCE):
		print "no bwa index file, make index."
		if not os.path.exists(BWA_INDEX):
			os.mkdir(BWA_INDEX)
		cmd_make_index=bwa_path+" index -p "+BWA_INDEX+'/'+GENOME + " -a bwtsw "+REFERENCE
		print cmd_make_index
		os.system(cmd_make_index)
	else:
		print "ERROR: no index file and no genome.fa to build it"
		os._exit(1)
	if os.path.exists(kallisto_index_path): 
		print "kallisto index path:%s"%(kallisto_index_path)
	else:
		print "ERROR: no index file and no genome.fa to build it"
		os._exit(1)
	#####check output directory###
	print "check output directory"
	if not os.path.exists(output_fold):
		os.mkdir(output_fold)
	if not os.path.exists(somatic_mutation_fold):
		os.mkdir(somatic_mutation_fold)
	if not os.path.exists(netmhc_out_fold):
		os.mkdir(netmhc_out_fold)
	if not os.path.exists(netctl_out_fold):
		os.mkdir(netctl_out_fold)
	if not os.path.exists(varscan_copynumber_fold):
		os.mkdir(varscan_copynumber_fold)
	if not os.path.exists(pyclone_fold):
		os.mkdir(pyclone_fold)
	if not os.path.exists(alignment_out_fold):
		os.mkdir(alignment_out_fold)
	if not os.path.exists(kallisto_out_fold):
		os.mkdir(kallisto_out_fold)
	if not os.path.exists(candidate_neoantigens_fold):
		os.mkdir(candidate_neoantigens_fold)
	if not os.path.exists(candidate_neoantigens_fold):
		os.mkdir(candidate_neoantigens_fold)			
	print "start fastq quality control"
	processes_0=[]
	q1=multiprocessing.Process(target=read_trimmomatic,args=(tumor_fastq_path_first,tumor_fastq_path_second,trimmomatic_path,adapter_path,tumor_fastq_clean_first,tumor_fastq_clean_second,tumor_fastq_unpaired_first,tumor_fastq_unpaired_second,))
	processes_0.append(q1)
	q2=multiprocessing.Process(target=read_trimmomatic,args=(normal_fastq_path_first,normal_fastq_path_second,trimmomatic_path,adapter_path,normal_fastq_clean_first,normal_fastq_clean_second,normal_fastq_unpaired_first,normal_fastq_unpaired_second,))
	processes_0.append(q2)
	for p in processes_0:
		p.daemon = True
		p.start()
	for p in processes_0:
		p.join()


	print "start stage 1"
	processes_1=[]
	if hla_str=="None":
		d1=multiprocessing.Process(target=hlatyping,args=(tumor_fastq_clean_first,tumor_fastq_clean_second,opitype_fold,opitype_out_fold,opitype_ext,prefix,))
 		#processes_1.append(d1)
 		#q.put('hlatyping') 		
 		#hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]
 	else:
 		print "hla type provided!"
 	#d2=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(normal_fastq_clean_first,normal_fastq_clean_second,'normal',CPU,BWA_INDEX,GENOME,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138_path,OneKG_path,mills_path,))
 	#processes_1.append(d2)
 	#q.put('normal_qc')
 	#d3=multiprocessing.Process(target=mapping_qc_gatk_preprocess,args=(tumor_fastq_path_first,tumor_fastq_path_second,'tumor',CPU,BWA_INDEX,GENOME,alignment_out_fold,prefix,REFERENCE,bwa_path,samtools_path,java_picard_path,GATK_path,dbsnp138_path,OneKG_path,mills_path,))
 	#processes_1.append(d3)
 	d4=multiprocessing.Process(target=kallisto_expression,args=(rna_fastq_1_path,rna_fastq_2_path,kallisto_path,kallisto_out_fold,prefix,kallisto_index_path,))
 	processes_1.append(d4)
 	#q.put('tumor_qc')
 	#for p in processes_1:
	#	p.daemon = True
	#	p.start()
	#for p in processes_1:
	#	p.join()
	#print 'stage 1 done.'
	if exp_file!="None" and os.path.exists(exp_file):
		print "check expression file done."
	elif exp_file=="None":
		print "no expression file provided."
	else:
		print "please check your expression file path!"
		os._exit(1)
	print 'start stage 2'
	processes_2=[]
	h1=multiprocessing.Process(target=varscan_somatic_caling_drift,args=(somatic_mutation_fold,alignment_out_fold,prefix,REFERENCE,vep_cache,samtools_path,varscan_path,vep_path,netmhc_out_fold,))
	processes_2.append(h1)
	#h2=multiprocessing.Process(target=indel_calling_drift,args=(strelka_out_fold,strelka_path,alignment_out_fold,prefix,REFERENCE,vep_cache,netmhc_out_fold,CPU,vep_path,))
	#processes_2.append(h2)
	#h3=multiprocessing.Process(target=varscan_copynumber_calling,args=(varscan_copynumber_fold,prefix,alignment_out_fold,REFERENCE,samtools_path,varscan_path,))
	#processes_2.append(h3)
	#for p in processes_2:
	#	p.daemon = True
	#	p.start()
	#for p in processes_2:
	#	p.join()
	if hla_str=="None":
 		hla_str=open(opitype_out_fold+'/'+prefix+"_optitype_hla_type").readlines()[0]
	print 'start stage 3'
	print exp_file
	processes_3=[]
	t2=multiprocessing.Process(target=varscan_neo,args=(snv_fasta_file,hla_str,snv_netmhc_out_file,netmhc_out_fold,split_num,prefix,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,))
	processes_3.append(t2)
	t3=multiprocessing.Process(target=indel_neo,args=(somatic_mutation_fold,prefix,vep_cache,netmhc_out_fold,vep_path,indel_fasta_file,hla_str,indel_netmhc_out_file,split_num,exp_file,binding_fc_aff_cutoff,binding_aff_cutoff,fpkm_cutoff,netctl_out_fold,netMHCpan_path,))
	processes_3.append(t3)
	#for p in processes_3:
	#	p.daemon = True
	#	p.start()
	#for p in processes_3:
	#	p.join()
	#print "stage 3 done."
	print 'start stage 4.'
	processes_4=[]
	l1=multiprocessing.Process(target=pyclone_annotation,args=(somatic_mutation_fold,varscan_copynumber_fold,prefix,pyclone_fold,netctl_out_fold,coverage,pyclone_path,cancer_type,))
	processes_4.append(l1)	
	#for p in processes_4:
	#	p.daemon = True
	#	p.start()
	#for p in processes_4:
	#	p.join()
	print 'stage 4 done.'
	print 'start stage 5.'
	processes_5=[]
	r1=multiprocessing.Process(target=immunogenicity_score_calculate,args=(final_neo_file,immunogenicity_score_ranking_file,))
	processes_5.append(r1)	
	r2=multiprocessing.Process(target=rankagg_score_calculate,args=(final_neo_file,neo_for_rankagg_file,neo_rank_agg_result,))
	processes_5.append(r2)
	#for p in processes_5:
	#	p.daemon = True
	#	p.start()
	#for p in processes_5:
	#	p.join()		
	print 'stage 5 done.'
	print 'Done!'
'''
	if os.path.exists(alignment_out_fold):
		shutil.rmtree(alignment_out_fold)
	if os.path.exists(somatic_mutation_fold):
		shutil.rmtree(somatic_mutation_fold)
	if os.path.exists(varscan_copynumber_fold):
		shutil.rmtree(varscan_copynumber_fold)
	if os.path.exists(netmhc_out_fold):
		shutil.rmtree(netmhc_out_fold)
	if os.path.exists(netctl_out_fold):
		shutil.rmtree(netctl_out_fold)
	if os.path.exists(kallisto_out_fold):
		shutil.rmtree(kallisto_out_fold)
	if os.path.exists(pyclone_fold):
		shutil.rmtree(pyclone_fold)
	if os.path.exists(opitype_out_fold):
		shutil.rmtree(opitype_out_fold)
	if os.path.exists(strelka_out_fold):
		shutil.rmtree(strelka_out_fold)
	shutil.rmtree(output_fold + '/' + 'STDOUT_summary.html')
'''
