# ======================================================================================================================
# Project: BioPipe
# Script : build_index.smk TODO check 
# Author : Peng Jia
# Date   : 2020.12.15
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
import os

configfile: "config.yaml"
##### config
dir_database = "/home/pengjia/reference_db/"
dir_log = "/home/pengjia/reference_db/logs/"

###### software
bwa = "/home/pengjia/miniconda3/envs/ngs/bin/bwa"
samtools = "/home/pengjia/miniconda3/envs/ngs/bin/samtools"
bowtie_build = "/home/pengjia/miniconda3/envs/tophat/bin/bowtie-build"
bowtie2_build = "/home/pengjia/miniconda3/envs/tophat/bin/bowtie2-build"
picard = "/home/pengjia/miniconda3/envs/ngs/bin/picard"
tophat2 = "/home/pengjia/miniconda3/envs/tophat/bin/tophat2"

#### control
##genome
dir_genome = dir_database + "genome/"
path_genome_input = config["reference"]["path"]
genome_prefix = config["species"] + "." + config["reference"]["version"]
path_genome_prefix = dir_genome + genome_prefix
path_genome = dir_genome + genome_prefix + ".fa"
path_genome_bwt_tag = path_genome + ".bwa.ok"
path_genome_bowtie_tag = path_genome + ".bowtie.ok"
path_genome_bowtie2_tag = path_genome + ".bowtie2.ok"
path_genome_picard = path_dict = path_genome.replace("fa", "dict").replace("fasta", "dict")

## gff and gtf
dir_annotation = dir_database + "annotation/"
path_gtf_input = config["annotation"]["gtf"]
path_gff_input = config["annotation"]["gff"]
path_gtf = dir_annotation + config["species"] + "." + config["annotation"]["version"] + ".gtf"
path_gff = dir_annotation + config["species"] + "." + config["annotation"]["version"] + ".gff"
path_genome_tophat_tag = path_genome + ".tophat2.ok"
path_genome_tophat_prefix = path_genome_prefix + ".tr"

aligner = ["bwa", "bowtie", "bowtie2", "star", "tophat2"]

path_genome_output = [path_dict]
if "tophat2" in aligner:
    aligner.append("bowtie")
    aligner.append("bowtie2")

if "bwa" in aligner:
    path_genome_output.append(path_genome_bwt_tag)

if "bowtie2" in aligner:
    path_genome_output.append(path_genome_bowtie2_tag)
if "bowtie" in aligner:
    path_genome_output.append(path_genome_bowtie_tag)
if "tophat2" in aligner:
    path_genome_output.append(path_genome_tophat_tag)

rule all:
    input:
         path_genome_output
# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO


rule load_genome:
    input:
         genome=path_genome_input,
    output:
          genome=path_genome,

    run:
        shell("cp {input.genome} {output.genome}")
        shell("{samtools} faidx {output}")

rule load_annotation:
    input:
         gtf=path_gtf_input,
         gff=path_gff_input,
    output:
          gtf=path_gtf,
          gff=path_gff,

    run:
        shell("cp {input.gtf} {output.gtf}")
        shell("cp {input.gff} {output.gff}")

rule BWA_INDEX:
    input:
         genome=path_genome,
    output:
          bwt_a=path_genome + ".amb",  # for bwa
          bwt_b=path_genome + ".bwt",  # for bwa
          bwt_c=path_genome + ".pac",  # for bwa
          bwt_d=path_genome + ".ann",  # for bwa
          bwt_e=path_genome + ".sa",  # for bwa
          tag=path_genome_bwt_tag
    threads: 1
    log:
       dir_log + "bwa_index/bwa_index.log"
    benchmark:
             dir_log + "bwa_index/bwa_index.tsv"
    run:
        shell("{bwa} index {input.genome}  2>{log} 1>{log}")
        shell("touch {output.tag}")

rule Bowtie2_Index:
    input:
         genome=path_genome,
    output:
          bwt_a=path_genome_prefix + ".1.bt2",
          bwt_b=path_genome_prefix + ".2.bt2",
          bwt_c=path_genome_prefix + ".3.bt2",
          bwt_d=path_genome_prefix + ".4.bt2",
          bwt_e=path_genome_prefix + ".rev.1.bt2",
          bwt_f=path_genome_prefix + ".rev.2.bt2",
          tag=path_genome_bowtie2_tag
    threads: 16
    log:
       dir_log + "bwa_index/bowtie2_index.log"
    benchmark:
             dir_log + "bwa_bowtie2/bowtie2_index.tsv"
    run:
        shell("{bowtie2_build} {input.genome} {path_genome_prefix} 2>{log} 1>{log}")
        shell("touch {output.tag}")

rule Bowtie_Index:
    input:
         genome=path_genome,
    output:
          bwt_a=path_genome_prefix + ".1.ebwt",
          bwt_b=path_genome_prefix + ".2.ebwt",
          bwt_c=path_genome_prefix + ".3.ebwt",
          bwt_d=path_genome_prefix + ".4.ebwt",
          bwt_e=path_genome_prefix + ".rev.1.ebwt",
          bwt_f=path_genome_prefix + ".rev.2.ebwt",
          tag=path_genome_bowtie_tag
    threads: 16
    log:
       dir_log + "bwa_index/bowtie2_index.log"
    benchmark:
             dir_log + "bowtie_index/bowtie2_index.tsv"
    run:
        shell("{bowtie_build} {input.genome} {path_genome_prefix} 2>{log} 1>{log}")
        shell("touch {output.tag}")

rule Picard_Index:
    input:
         path_genome
    output:
          path_dict  # for gatk
    log:
       dir_log + "picard_index/picard_index.log"
    benchmark:
             dir_log + "picard_index/picard_index.tsv"
    run:
        shell("{picard} CreateSequenceDictionary R={input} O={output}  2>{log} 1>{log} ")

rule Tophat2_Index:
    input:
         fasta=path_genome,
         gtf=path_gtf
    output:
          bwt_a=path_genome_tophat_prefix + ".1.bt2",
          tag=path_genome_tophat_tag
    log:
       dir_log + "tophat2_index/tophat2_index.log"
    threads: 16
    benchmark:
             dir_log + "tophat2_index/tophat2_index.tsv"
    run:
        fasta = str(input.fasta).rstrip(".fa").rstrip(".fasta")
        shell("{tophat2} -G {input.gtf} -p {threads} --transcriptome-index={path_genome_prefix}.tr {fasta} "
              " 2>{log} 1>{log}")
        shell("touch {output.tag}")
# shell("touch {output.tag}")
# shell("{picard} CreateSequenceDictionary R={input} O={output}  2>{log} 1>{log} ")
