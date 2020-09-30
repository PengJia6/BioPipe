# ======================================================================================================================
# Project: DNAseq
# Script : var_manta.smk
# Author : Peng Jia
# Date   : 2020.07.19
# Email  : pengjia@stu.xjtu.edu.cn
# Description: SV calling with manta
# ======================================================================================================================

# localrules: a
# ruleorder: a > b

# ======================================================================================================================
# rules: Manta_conf
# description: Manta configure
# input: bam and fasta
# output: workflow
rule Manta_conf:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
    output:
          prefix=directory(path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}_manta"),
          tag=temp(path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}"),
          # workflow=path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}_manta/runWorkflow.py",
          # vcfgz=path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}/{bam_sample}-smoove.genotyped.vcf.gz",
          # final_vcfgz=path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}.smoove.raw.vcf.gz"
    log:
       path_log + "gremlineVar/manta/perSample/{bam_sample}/{bam_sample}.manta_conf.logs"
    benchmark:
             path_bm + "gremlineVar/manta/perSample/{bam_sample}/{bam_sample}.manta_conf.tsv"

    threads: config["threads"]["Manta_conf"]
    params:
          extra="",
    run:
        shell("{path_manta}configManta.py --bam {input.bam} --referenceFasta {input.ref} "
              "--runDir {output.prefix} 2>{log} 1>{log}")
        shell("touch {output.tag}")
# shell("command 2")
# ======================================================================================================================
# rules: Manta
# description:
# input: configure workflow
# output: vcf.gz
rule Manta:
    input:
        prefix=path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}",
    output:
          tmp_vcf=path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}_manta/results/variants/diploidSV.vcf.gz",
          vcfgz=path_data + "germlineVar/manta/perSample/{bam_sample}/{bam_sample}_manta.raw.vcf.gz",

    log:
       path_log + "gremlineVar/manta/perSample/{bam_sample}/{bam_sample}.manta.logs"
    benchmark:
             path_bm + "gremlineVar/manta/perSample/{bam_sample}/{bam_sample}.manta.tsv"

    threads: config["threads"]["Manta"]
    params:
          extra="",
    shell:
         """
         export PATH={path_manta}:$PATH
         chmod +x {input.prefix}_manta/runWorkflow.py
         {input.prefix}_manta/runWorkflow.py -j {threads} 2>{log} 1>{log}
         ln {output.tmp_vcf} {output.vcfgz}
         ln {output.tmp_vcf}.tbi {output.vcfgz}.tbi
         
         """
