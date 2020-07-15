# ==============================================================================
# Project: DNAseq
# Script : var_lumpy.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Variants calling with lumpy
# ==============================================================================

# localrules: a
# ruleorder: a > b


# rules: Lumpy_Get_Discordants:
# description: get Discordants reads
# input: bam
# output: bam
rule Lumpy_Get_Discordants:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
         excl=path_genome_prefix+"exclude.cnvnator.bed"
    output:
        bam=temp(path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}.discordants.bam"),
    log:
       ""
    benchmark:
             ""
    threads: 8
    params:
          extra="",
    run:
        shell("{path_samtools}samtools view -@ {threads} -b -F 1294 {input.bam}| "
               "{path_samtools}samtools sort  -@ {threads} -o {output.bam} 2>{log} 1>{log}")

# rules: TODO
# description: TODO
# input: TODO
# output: TODO
rule TODO2:
    input:
         "",
    output:
          "",
    log:
       ""
    benchmark:
             ""
    threads: 8
    params:
          extra="",
    shell:
         """
         command 1 
         command 2
         """