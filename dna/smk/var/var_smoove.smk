# ======================================================================================================================
# Project: DNAseq
# Script : var_smoove.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Variant Calling by Lumpy
# Homepage smoove: https://github.com/brentp/smoove
# Homepage lumpy: https://github.com/arq5x/lumpy-sv
# ======================================================================================================================


# ======================================================================================================================
# rules: Smoove
# description: Variant Calling using Smoove (best practice of Lumpy)
# input: bam
# output: vcf
rule Smoove:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
         excl=path_genome_prefix+"exclude.cnvnator.bed"
    output:
          prefix=directory(path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}"),
          # vcfgz=path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}/{bam_sample}-smoove.genotyped.vcf.gz",
          final_vcfgz=path_data + "germlineVar/smoove/perSample/{bam_sample}/{bam_sample}.smoove.raw.vcf.gz"

    log:
       path_log + "gremlineVar/smoove/perSample/{bam_sample}/{bam_sample}.smoove.logs"
    benchmark:
             path_bm + "gremlineVar/smoove/perSample/{bam_sample}/{bam_sample}.smoove.tsv"
    threads: config["threads"]["Smoove"]
    params:
          extra="",
    shell:
        """
        export PATH={path_lumpy}:$PATH
        {path_smoove}smoove call -x -d --genotype --name {wildcards.bam_sample} --outdir {output.prefix} --exclude {input.excl} --fasta {input.ref} -p {threads} {input.bam} 2>{log} 1>{log}
        mv {output.prefix}/{wildcards.bam_sample}-smoove.genotyped.vcf.gz {output.final_vcfgz}
        """
# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
# rule TODO2:
#     input:
#          "",
#     output:
#           "",
#     log:
#        ""
#     benchmark:
#              ""
#     threads: 8
#     params:
#           extra="",
#     shell:
#          """
#          command 1
#          command 2
#          """