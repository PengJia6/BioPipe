# ======================================================================================================================
# Project: DNAseq
# Script : var_delly.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Variants detection by delly
# Homepage: https://github.com/dellytools/delly
# ======================================================================================================================
# rules: Delly_Call
# description: SV Calling
# input: sorted bam
# output: bcf
rule Delly_Call:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
         excl=path_genome_prefix+"exclude.cnvnator.bed"
    output:
          bcf=path_data + "germlineVar/delly/perSample/{bam_sample}/{bam_sample}.delly.call.bcf",

    log:
       path_log + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellycall.logs"
    benchmark:
             path_bm + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellycall.tsv"
    threads: config["threads"]["Delly_Call"]
    params:
          extra="",
    run:
        shell("{path_delly}delly call -g {input.ref} -x {input.excl} -o {output.bcf} {input.bam} 2>{log} 1>{log}")
# ======================================================================================================================
# rules: Delly_GT
# description: Genotype the delly call result
# input: bcf and bam
# output: genotyped vcf
rule Delly_GT:
    input:
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai",
         bcf=rules.Delly_Call.output.bcf,
         excl=path_genome_prefix+"exclude.cnvnator.bed"

    output:
          bcf=path_data + "germlineVar/delly/perSample/{bam_sample}/{bam_sample}.delly.raw.bcf",
    log:
       path_log + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellyGT.logs"

    benchmark:
             path_bm + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellyGT.tsv"

    threads: config["threads"]["Delly_GT"]
    params:
          extra="",
    run:
        shell("{path_delly}delly call -g {input.ref} -x {input.excl} -o {output.bcf} -v {input.bcf} {input.bam} 2>{log} 1>{log}")
# ======================================================================================================================
# rules: Delly2Bcf2vcfgz
# description: covert bcf file to vcf.gz
# input: bcf
# output: vcf.gz
rule Delly2Bcf2vcfgz:
    input:
         bcf=rules.Delly_GT.output.bcf
    output:
          vcfgz=path_data + "germlineVar/delly/perSample/{bam_sample}/{bam_sample}.delly.raw.vcf.gz",
    log:
       path_log + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellybcf2vcfgz.logs"

    benchmark:
             path_bm + "gremlineVar/delly/perSample/{bam_sample}/{bam_sample}.dellybcf2vcfgz.tsv"

    threads: config["threads"]["Delly2Bcf2vcfgz"]
    params:
          extra="",
    run:
        shell("{path_bcftools}bcftools view {params.extra} -Oz -o {output.vcfgz} {input} 2>{log} 1>{log}")
