# Next generation sequencing analysis pipeline 
## How to use the pipeline
  1 copy the directory to you path 
      ```
       cp -r /path/to/BipPipe/ngs_germline /mypath/BioPipe
       cd /mypath/BioPipe/ngs_germline

      ```
      
   or clone the pipeline to your path  
      ```
       git clone xxxx   
       BioPipe/ngs_germline   
      ```
  2. configure your envelopment
   * config the input data in conf/sample_info.csv  
   such as 

   |sample |unit|  LB  |  PL |  fq1  |fq2  |  
   | ---- |  ----|---- |  ----|---- |  ----|   
   | HG001 | L1|L1 | ILM|HG001_L1.R1.fq | HG001_L1.R2.fq|
   | HG001 | L2|L2 | ILM|HG001_L2.R1.fq | HG001_L2.R2.fq|
   | HG003 | L1|L1 | MGI|HG003_L1.R1.fq | HG003_L1.R2.fq|
   
   sample: sample name 
   unit: library unit 
   LBL: library name 
   PL: sequencing platform, ILM/BGI/MGI
   fq1: path of the read1 for paired-end sequencing 
   fq2: path of the read2 for paired-end sequencing 
   * if you have high quality
    sequencing config the input data in conf/sample_info.csv
      
   |case|path|  
   | ---- |  ----|   
   | HG001 | /path/to/HG001.HQ.bam|
   | HG003 | /path/to/HG003.HQ.bam|
  

     
   
  
     







# DNAseq
## Done
* add passqc (√)
* make test samples (√)
* add pigz rules for gzip (√)
* Process error and output of pbs（√）  
* remove the dependence of wrapper, use run and shell （√）
* add more caller for pipeline （√）

## Doing

* add qc report 
* add bam report 
---
### QC 
* fastqc  (√)
* fastp   (√)
* trim..  
* multiqc 

### aligner

* bwa (√)
* bowtie 
* bowtie2 

### caller 
#### germline 
* gatk hc (√) 
* samtools （√）
* freebayse （√）
* varscan （√）
* ...

#### somatic 
