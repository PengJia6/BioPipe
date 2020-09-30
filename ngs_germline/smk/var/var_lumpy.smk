# ==============================================================================
# Project: DNAseq
# Script : var_lumpy.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Variants calling with lumpy
# Homepage lumpy: https://github.com/arq5x/lumpy-sv
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
    output:
          bam=path_data + "germlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.discordants.bam",
    log:
       path_log + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.discordants.logs"
    benchmark:
             path_bm + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.discordants.tsv"

    threads: config["threads"]["Lumpy_Get_Discordants"]
    params:
          extra="",
    run:
        shell("{path_samtools}samtools view -@ {threads} -b -F 1294 {input.bam}| "
              "{path_samtools}samtools sort  -@ {threads} -o {output.bam} 2>{log} 1>{log}")

# rules: Lumpy_Get_SplitReads
# description: get Discordants reads
# input: bam file
# output: splitters reads
rule Lumpy_Get_SplitReads:
    input:
         unpack(getHQbamsample),
    output:
          bam=path_data + "germlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.splitters.bam",
    log:
       path_log + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.splitters.logs"
    benchmark:
             path_bm + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.splitters.tsv"

    threads: config["threads"]["Lumpy_Get_SplitReads"]
    params:
          extra="",
    run:
        shell("{path_samtools}samtools view -h -@ {threads} {input.bam}| "
              "{path_lumpy}extractSplitReads_BwaMem -i stdin | "
              "{path_samtools}samtools view -@ {threads} -Sb | "
              "{path_samtools}samtools sort  -@ {threads} -o {output.bam} 2>{log} 1>{log}")
# rules: Lumpyexpress
# description: lumpy calling
# input: bam file
# output: splitters reads
rule Lumpyexpress:
    input:
         unpack(getHQbamsample),
         discordants=rules.Lumpy_Get_Discordants.output.bam,
         split=rules.Lumpy_Get_SplitReads.output.bam,
         ref=path_genome,
         sindex=path_genome + ".fai",
         excl=path_genome_prefix + "exclude.cnvnator.bed"
    output:
          vcf=path_data + "germlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.lumpy.raw.vcf",
          tmp=directory(path_data + "germlineVar/lumpy/perSample/{bam_sample}/{bam_sample}_tmp"),
          vcfgz=path_data + "germlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.lumpy.raw.vcf.gz",

    log:
       path_log + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.lumpy.logs"
    benchmark:
             path_bm + "gremlineVar/lumpy/perSample/{bam_sample}/{bam_sample}.lumpy.tsv"

    threads: config["threads"]["Lumpyexpress"]
    params:
          extra="",
    shell:
         """
         export PATH={path_lumpy}:$PATH
         mkdir {output.tmp}
         {path_lumpy}lumpyexpress -R {input.ref} -B {input.bam} -S {input.split} -D {input.discordants} -o {output.vcf} -x {input.excl} -T {output.tmp} 2>{log} 1>{log}
           {path_bgzip}bgzip -c {output.vcf} > {output.vcfgz} 
         """
