#!/usr/bin/python
# -*- coding: UTF-8 -*-
###########netMHC result parsing and filter based on binding affinity and FPKM #########
import sys,getopt,os
import pandas as pd 


opts,args=getopt.getopt(sys.argv[1:],"hi:g:e:o:b:a:f:s:",["input_netmhc_file","input_fasta","expression_fpkm_file","out_dir","binding_affinity_cutoff","binding_affinity_foldchange_cutoff","fpkm_cutoff","sample_id"])
input_netmhc_file=""
input_fasta=""
expression_fpkm_file=""
out_dir=""
binding_affinity_cutoff=500
binding_affinity_foldchange_cutoff=1
fpkm_cutoff=1
sample_id=""
USAGE='''usage: python netMHC_result_parse.py -i <input_netmhc_file> -g <input_fasta> -o <outdir> -s <sample_id> [option]
		required argument:
			-i | --input_netmhc_file : input file,result from netMHC
			-g | --input_fasta : input fasta file for netMHC
			-o | --out_dir : output directory
			-s | --sample_id : sample id
		optional argument:
			-e | --expression_fpkm_file : expression profile with fpkm value
			-b | --binding_affinity_cutoff : pipetide binding affinity cutoff , default: 500
			-a | --binding_affinity_foldchange_cutoff : pipetide binding affinity fold change between mutate peptide and wild type peptide cutoff , default: 1
			-f | --fpkm_cutoff : FPKM value cutoff, default :1'''
	
for opt,value in opts:
	if opt =="h":
		print USAGE
		sys.exit(2)
	elif opt in ("-i","--input_netmhc_file"):
		input_netmhc_file=value
	elif opt in ("-g","--input_fasta"):
		input_fasta=value
	elif opt in ("-e","--expression_fpkm_file"):
		expression_fpkm_file =value
	elif opt in ("-o","--out_dir"):
		out_dir =value  
	elif opt in ("-b","--binding_affinity_cutoff"):
		binding_affinity_cutoff =value
	elif opt in ("-a","--binding_affinity_foldchange_cutoff"):
		binding_affinity_foldchange_cutoff =value
	elif opt in ("-f","--fpkm_cutoff"):
		fpkm_cutoff =value 
	elif opt in ("-s","sample_id"):
		sample_id=value
#print coverage
if (input_netmhc_file =="" or input_fasta =="" or out_dir =="" or sample_id==""):
	print USAGE
	sys.exit(2)		
#######extract full animo acid change##
Full_header=[]
Full_gene=[]
with open(input_fasta) as f:
	data=f.read()
for line in data.strip().split('\n'):
	if line.startswith('>WT'):
		full_aa_change=line.split('_')[-1]
		full_gene=line.split('_')[1]
		Full_header.append(full_aa_change)
		Full_gene.append(full_gene)
	else:
		continue
dup_full_header=[]
dup_full_gene=[]
hla_num=6
i=0
while i<6:
	for j in range(len(Full_header)):
		dup_full_header.append(Full_header[j])
		dup_full_gene.append(Full_gene[j])
	i=i+1


#print Full_header
#print len(dup_full_header)
#print dup_full_header[2625]
######## extract candidate neoantigens####
#print expression_fpkm_file
with open(input_netmhc_file) as f:
    data = f.read()

nw_data = data.split('-----------------------------------------------------------------------------------\n')
WT_header = []
MT_header = []
WT_neo = []
MT_neo = []
for i in range(len(nw_data)):
    if i%8 == 3:
        wt_pro_name = nw_data[i].strip('\n').split('.')[0]
        WT_header.append(wt_pro_name)
    elif i%8 == 2:
        wt_neo_data = nw_data[i].strip().split('\n')
        WT_neo.append(wt_neo_data)
    elif i%8 == 7:
        mt_pro_name = nw_data[i].strip('\n').split('.')[0]
        MT_header.append(mt_pro_name)
    elif i%8 == 6:
        mt_neo_data = nw_data[i].strip().split('\n')
        MT_neo.append(mt_neo_data)
WB_SB_MT_record = []
WT_record = []
aa_record=[]
gene_record=[]
#print MT_neo

for i in range(len(MT_neo)):
	for j in range(len(MT_neo[i])):
		if MT_neo[i][j].endswith('WB') or MT_neo[i][j].endswith('SB'):
			#print i,j
			aa_record.append(dup_full_header[i])
			gene_record.append(dup_full_gene[i])
			WB_SB_MT_record.append(MT_neo[i][j])
			WT_record.append(WT_neo[i][j])


candidate_neo = open(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt",'w')
candidate_neo.write('\t'.join(['#HLA_type','Gene','AA_change','MT_pep','WT_pep','MT_Binding_affinity','WT_Binding_affinity','MT_Binding_level','WT_Binding_level','MT_Binding_level_des','WT_Binding_level_des','fold_change','DAI']) + '\n')
for i in range(len(WB_SB_MT_record)):
    mt_record = [line for line in WB_SB_MT_record[i].split(' ') if line!='']
    HLA_tp = mt_record[1]
    gene = gene_record[i]
    ani_change = aa_record[i]
    mt_pep = mt_record[2]
    mt_binding_aff = mt_record[12]
    mt_binding_level=mt_record[13]
    mt_binding_level_des = mt_record[-1]
    wt_record = [i for i in WT_record[i].split(' ') if i!='']
    wt_pep = wt_record[2]
    wt_binding_aff = wt_record[12]
    wt_binding_level=wt_record[13]
    if wt_record[-1]=='SB' or wt_record[-1]=='WB':
        wt_binding_level_des = wt_record[-1]
    else:
        wt_binding_level_des = 'NB'
    fold_change = float(mt_binding_aff)/float(wt_binding_aff)
    DAI = float(mt_binding_aff) - float(wt_binding_aff)
    out_line = '\t'.join((HLA_tp,gene,ani_change,mt_pep,wt_pep,mt_binding_aff,wt_binding_aff,mt_binding_level,wt_binding_level,mt_binding_level_des,wt_binding_level_des,str(fold_change),str(DAI)))
    candidate_neo.write(out_line + '\n')
candidate_neo.close()
    



######neoantigens filtering#####
##including Binding affinity ,localized peptide, multiple length peptide screen, differential AI> ##, neoantigens stability, gene FPKM>1 #####
data=pd.read_table(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt",header=0,sep='\t')
if expression_fpkm_file=='no_exp':
	final_filter_data=data[(data.MT_Binding_affinity<int(binding_affinity_cutoff)) & (data.fold_change<int(binding_affinity_foldchange_cutoff))] 
elif os.path.exists(expression_fpkm_file):
	first_filter_data=data[(data.MT_Binding_affinity<int(binding_affinity_cutoff)) & (data.fold_change<int(binding_affinity_foldchange_cutoff))]
	exp = pd.read_table(expression_fpkm_file,header=0,sep='\t')
	gene_exp = exp.loc[:,['Gene Name','TPM']]
	neo_merge_exp = pd.merge(first_filter_data,gene_exp,left_on='Gene',right_on='Gene Name',how='left')
	#print neo_merge_exp
	final_filter_data=neo_merge_exp#[neo_merge_exp.TPM>float(fpkm_cutoff)]
else:
	print "could not find expresion file,check if the file exists!"
	sys.exit(2)
#os.remove(out_dir+'/'+sample_id+"_tmp_neo_candidate.txt")
print final_filter_data
final_filter_data.to_csv(out_dir+'/'+sample_id+"_final_neo_candidate.txt",header=1,sep='\t',index=0)

    







            
