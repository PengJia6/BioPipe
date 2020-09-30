# ==============================================================================
# Project: DNAseq
# Script : var_breakdancer.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Breakdancer Calling
# https://github.com/genome/breakdancer
# ==============================================================================

# localrules: a
# ruleorder: a > b


# rules: BreakDancer_Conf
# description: Build config file for breakdancer
# input: bam
# output: config
rule BreakDancer_Conf:
    input:
         unpack(getHQbamsample),
    output:
          conf=path_data + "germlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.breakdancer.conf",
    log:
       path_log + "gremlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.breakdancer_conf.logs"
    benchmark:
             path_bm + "gremlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.breakdancer_conf.tsv"

    threads: config["threads"]["BreakDancer_Conf"]
    params:
          extra="",
    shell:
        """
        export PATH={path_breakdancer}:$PATH
        {path_breakdancer}bam2cfg.pl {input.bam} > {output.conf}
        """

# rules: BreakDancer
# description: Breakdancer calling
# input: conf
# output: breakdancer out
rule BreakDancer:
    input:
         rules.BreakDancer_Conf.output.conf
    output:
          conf=path_data + "germlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.raw.breakdancer",
    log:
       path_log + "gremlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.breakdancer.logs"
    benchmark:
             path_bm + "gremlineVar/breakdancer/perSample/{bam_sample}/{bam_sample}.breakdancer_conf.tsv"

    threads: config["threads"]["BreakDancer"]
    params:
          extra="",
    shell:
         """
         export PATH={path_breakdancer}:$PATH
         {path_breakdancer}breakdancer-max {input} > {output}
         """
