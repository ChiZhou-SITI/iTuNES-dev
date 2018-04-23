# iTuNES User Manual

Chi Zhou  
Biological and Medical Big data Mining Lab  
Tongji University  
Shanghai, China. 

## Table of Contents
1. [General Description](#general-description)  
2. [Dependencies](#dependencies)  
    - [Required software](#required-software)  
    - [Python packages](#python-packages)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [Input Files](#input-files)  
    - [VCF file](#vcf-file)  
    - [Expression file](#expression-file)  
    - [References](#references)  
6. [Output Files](#output-files)  
    - [Column explanation](#column-explanation)  
7. [Test Example](#test-example)  
    - [Data preparation](#data-preparation)  
        * [Recommended preprocessing of next generation sequencing (NGS)
data](#recommended-preprocessing-of-next-generation-sequencing-(ngs)-data)  
        * [Data cleanup](#data-cleanup)
        * [WXS data](#wxs-data)
        * [RNAseq](#rnaseq)
        * [HLA typing](#hla-typing)  

## General Description

Given matched tumor-normal whole exome sequencing and tumor RNA-seq sequencing data as input, iTuNes infers HLA sub-types, mutated peptides (neo-peptide), variant allele frequency, expression profile etc feature information. Based on these feature, a model- based refined ranking-score scheme could identify which of the neo-peptides have strong immunogecity.

## Dependencies  


#### Hardware:
iTuNEs currently test on x86_64 on ubuntu 16.04.

#### Required software:
* [Python 2.7](https://www.python.org/downloads/release/python-2712/)
* [R 3.2.3](https://cran.r-project.org/src/base/R-3/R-3.2.3.tar.gz)
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](https://github.com/Ensembl/ensembl-vep)
* [bwa](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools)
* [strelka](https://github.com/Illumina/strelka)
* [opitype](https://github.com/FRED-2/OptiType)
* [pyclone](https://bitbucket.org/aroth85/pyclone/wiki/Tutorial)
* [GATK 3.7](https://software.broadinstitute.org/gatk/best-practices/)
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)
* [Varscan2](http://varscan.sourceforge.net/)
* [kallisto](http://pachterlab.github.io/kallisto/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [vcftools](http://vcftools.sourceforge.net/)
* [blast](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)


#### Python modules:
    multiprocessing
    pyper
    sklearn

## Installation 

1. Install all software listed above. 
2. Install multiprocessing, pyper and sklearn with the following command:

        pip install multiprocessing
        pip install pyper
        pip install sklearn


3. Download or clone the iTuNEs repository to your local system:

        git clone https://github.com/XIAOCHIZI/iTuNES-dev.git

4. Reference data includes genome fasta, cDNA, peptide, cosmic reference(GRCh38 build) could be downloaded through:

        bash data_download.sh
        
    all reference data would be in the fold `database`, including:

        [Genome reference]
        Homo_sapiens_assembly38.fasta
        dbsnp_138.hg38.vcf.gz
        hapmap_3.3.hg38.vcf.gz
        1000G_phase1.snps.high_confidence.hg38.vcf.gz
        1000G_omni2.5.hg38.vcf.gz
        Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        CosmicCodingMuts_chr_M_sorted.vcf.gz
        
        [EnsemblVEP]
        homo_sapiens_vep_89_GRCh38.tar.gz
        
        [cDNA and protein]
        Homo_sapiens.GRCh38.cdna.all.fa
        Homo_sapiens.GRCh38.pep.all.fa

5. Fill in the `config.yaml` file with your local path, make sure you have installed all above software and have downloaded
reference data.You should be aware that the version of VEP library you use should match the references used (peptide and cDNA). E.g. in the example above is used version/release 85 of GRCh38.


## Usage

After installation iTuNES is called as follows.The user should  
The config file is specified using the `-c` option

    path/to/iTuNES.py -c path/to/config.yaml

a detailed explaination is in example `config.yaml` file, you should replace the path of required software and reference file in your system.

## Input Files 

MuPeXI accepts a VCF file of somatic mutation calls optimally obtained from either MuTect
or MuTect2. The VCF files do not need to be modified; the "raw" output VCF file can
be put directly into MuPeXI. 

### VCF file 

Compact example of a VCF file:

        ##fileformat=VCFv4.2
        ##GATKCommandLine.MuTect2=<ID=MuTect2,Version=3.5-0-g36282e4
        ##SAMPLE=<ID=NORMAL,SampleName=TCGA-XV-A9W5_N
        ##SAMPLE=<ID=TUMOR,SampleName=TCGA-XV-A9W5_T
        ##reference=file:///home/projects/pr_46630/data/references/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fa
        #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
        chr1    948711  .   C   G   .   germline_risk . . . .
        chr1    1657358 .   T   TA  .   alt_allele_in_normal . . . .
        chr1    1986752 rs4233028   A   G   .   germline_risk . . . .   
        chr1    3431704 rs2493274   G   C   .   germline_risk . . . .  
        chr1    3631978 rs2244942   T   C   .   germline_risk . . . .  
        chr1    3839305 rs1891940   T   C   .   clustered_events;germline_risk  . . . .

A full example of a VCF file can be found on the MuPeXI webserver
[here](http://www.cbs.dtu.dk/services/MuPeXI/example.vcf). 

### Expression file

It is optional, but preferable, to provide a file with expression values as input to add
the expression of each transcript where a mutated peptide was extracted.
The expression files used for testing MuPeXI were generated from raw RNA-seq data using
Kallisto. The files should be tab separated and include Ensembl transcript ID (ENST) and
mean expression WITHOUT a header. 

        ENST00000456328.2   0.868567715
        ENST00000450305.2   0
        ENST00000488147.1   2.72373575
        ENST00000619216.1   0
        ENST00000473358.1   0
        ENST00000469289.1   0
        ENST00000607096.1   0

A full example of an expression file can be found on the MuPeXI webserver
[here](http://www.cbs.dtu.dk/services/MuPeXI/example_expression.tsv).

It should be noted that MuPeXI takes both expression values determined on transcript and
gene level, though transcript is preferable. If gene level is used (ENSG...) the `-E gene`
option should be used. 

### References 
The following references are required for MuPeXI to run:
* Peptide  
    The peptide reference is a FASTA file containing all peptides of the human proteome.
* cDNA  
    The cDNA reference is a FASTA file containing all cDNA sequences of the human
    proteome.  

These references can be acquired from the
[Ensembl website](http://www.ensembl.org/Homo_sapiens/Info/Index).  
The most recent release is found under Gene annotation > Download genes, cDNAs, ncRNA,
protein (FASTA) > pep (for peptide reference) and > cdna (for cDNA reference).  
It should be emphasized that it is of very high importance that the references and VEP
match in release version (e.g. release-85).


The following reference are optional but preferable:
* Cosmic  
        TSV file containing known cancer driver genes. The cancer gene census can be
        downloaded from the [COSMIC](http://cancer.sanger.ac.uk/census) website.  

