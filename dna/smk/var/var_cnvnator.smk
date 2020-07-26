# ==============================================================================
# Project: DNAseq
# Script : var_cnvnator.smk TODO check 
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# https://github.com/abyzovlab/CNVnator
# ==============================================================================

# localrules: a
# ruleorder: a > b


# rules: CNVnator
# description: TODO
# input: TODO
# output: TODO
rule CNVnator:
    input: 
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
         root=path_data + "germlineVar/cnvnator/perSample/{bam_sample}/{bam_sample}.cnvnator.root",
         out=path_data + "germlineVar/cnvnator/perSample/{bam_sample}/{bam_sample}.raw.cnvnator"

    log:
       path_log + "gremlineVar/cnvnator/perSample/{bam_sample}/{bam_sample}.cnvnator.logs"
    benchmark:
             path_bm + "gremlineVar/cnvnator/perSample/{bam_sample}/{bam_sample}.cnvnator.tsv"
    threads: config["threads"]["CNVnator"]
    params:
          extra="",
    shell:
        """
        {path_cnvnator}cnvnator -root {output.root} -tree {input.bam} -unique 2>{log} 1>{log}
        {path_cnvnator}cnvnator -root {output.root} -d {input.ref} -his 100 2>{log} 1>{log}
        {path_cnvnator}cnvnator -root {output.root} -stat 100 2>{log} 1>{log}
        {path_cnvnator}cnvnator -root {output.root} -partition 100 2>{log} 1>{log}
        {path_cnvnator}cnvnator -root {output.root} -call 100 > {output.out}
        """
        # shell("command 2")

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