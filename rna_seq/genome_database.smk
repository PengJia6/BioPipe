# ======================================================================================================================
# Project: BioPipe
# Script : genome_database.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

# localrules: a
# ruleorder: a > b
include: "conf/input_and_output.smk"
# ======================================================================================================================
# rules: STAR_INDEX
# description: Generate genome index for STAR
# input: genome
# output: INDEX
rule STAR_Index:
    input:
         fasta=path_genome,
         gtf=path_gtf,
         gff=path_gff,
    output:
          index_dir=directory(path_STAR_index + "_{read_len}"),
          index_tag=path_STAR_index + "_{read_len}_tag"
    log:
       path_log + "genome/star_index_{read_len}.log"
    benchmark:
             path_log + "genome/star_index_{read_len}.tsv"
    threads: config["threads"]["STAR_Index"]
    params:
          extra="",
    run:
        shell(path_star + " --runThreadN {threads} --runMode genomeGenerate "
                          "--genomeDir {output.index_dir} --genomeFastaFiles {input.fasta} "
                          "--sjdbGTFfile {input.gtf} --sjdbOverhang {wildcards.read_len} 2>{log} 1>{log}")
        shell("touch {output.index_tag}")
# shell("command 2")

# ======================================================================================================================

rule tophat2_index:
    input:
         fasta=path_genome,
         gtf=path_gtf,
         gff=path_gff,
    output:
          tag=path_tophat2_index,
    log:
       path_log + "genome/tophat2_index.log"
    benchmark:
             path_log + "genome/tophat2_index.tsv"
    threads: config["threads"]["Tophat2_Index"]
    params:
          extra="",
    run:
        shell("ln -sf {input.fasta} {output.tag}.fa")
        shell("{path_tophat2} -G {input.gtf} --transcriptome-index={output.tag} {output.tag} "
              " 2>{log} 1>{log}")
        shell("touch {output.tag}")
