# ======================================================================================================================
# Project: DNAseq
# Script : genome.smk TODO check 
# Author : Peng Jia
# Date   : 2020.07.21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

# localrules: a
# ruleorder: a > b

# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
path_genome = str(path_data + "genome/" + config["ref"]["name"] + "/" + config["ref"]["genome"].split("/")[-1])
path_dict = path_genome.replace("fa", "dict").replace("fasta", "dict")
path_dict_orginal = config["ref"]["genome"].replace("fa", "dict").replace("fasta", "dict")
path_genome_list = [path_genome + ".fai", path_genome + ".bwt", path_dict, path_dict + "_chrom"]
path_orginal_genome_prefix = config["ref"]["genome"].rstrip("fa").rstrip("fasta")
path_genome_prefix = str(path_data + "genome/" + config["ref"]["name"] + "/" +
                         config["ref"]["genome"].split("/")[-1].rstrip("fa").rstrip("fasta"))
rule LoadGenome:
    input:
         config["ref"]["genome"]
    output:
          path_genome
    shell:
         """
         ln -sr {input} {output}
         touch -h {output}
         """
rule LoadGenomeFile:
    input:
        sv_excul=path_orginal_genome_prefix+"exclude.cnvnator.bed"
    output:
        sv_excul=path_genome_prefix+"exclude.cnvnator.bed"
    shell:
         """
         ln -sr {input.sv_excul} {output.sv_excul}
         touch -h {output.sv_excul}
         """



rule GenomeIndexSamtools:
    input:
         path_genome
    output:
          fai=path_genome + ".fai",  # samtools faidx
          fai_chrom=path_genome + ".fai_chrom",  # samtools faidx

    log:
       path_log + "genomeindex/genomeindex_samtools.log"
    run:
        if os.path.exists(config["ref"]["genome"] + ".fai"):
            shell("ln -sr {orign_ref}.fai {path_genome}.fai".format(orign_ref=config["ref"]["genome"],
                                                                    path_genome=path_genome))
        else:
            shell("{path_samtools}samtools faidx {input}  2>{log} 1>{log} ")

        if os.path.exists(config["ref"]["genome"] + ".fai_chrom"):
            shell("ln -sr {orign_ref}.fai_chrom {path_genome}.fai_chrom".format(orign_ref=config["ref"]["genome"],
                                                                                path_genome=path_genome))
        else:
            shell("head -n {num} {path_genome}.fai > {path_genome}.fai_chrom".format(num=config["ref"]["chrom_num"],
                                                                                     path_genome=path_genome))

rule GenomeIndexBwa:
    input:
         path_genome
    output:
          path_genome + ".amb",  # for bwa
          path_genome + ".bwt",  # for bwa
          path_genome + ".pac",  # for bwa
          path_genome + ".ann",  # for bwa
          path_genome + ".sa",  # for bwa
    log:
       path_log + "genomeindex/genomeindex_bwa.log"
    run:
        bwt_suffix = ["amb", "ann", "bwt", "pac", "sa"]
        bwt_stat = []
        for i in bwt_suffix:
            bwt_stat.append(os.path.exists(config["ref"]["genome"] + "." + i))
        if (False in bwt_stat):
            shell("{path_bwa}bwa index {input}  2>{log} 1>{log}")
        else:
            for i in bwt_suffix:
                shell("ln -sr " + config["ref"]["genome"] + "." + i + " " + path_genome + "." + i)
                shell("touch -h {path_genome}.{i}")

# TODO check
rule GenomeIndexPicard:
    input:
         path_genome
    output:
          path_dict  # for gatk
    log:
       path_log + "genomeindex/genomeindex_picard.log"
    run:
        if os.path.exists(path_dict_orginal):
            shell("ln -sr {path_dict_orginal} {path_dict}")
        else:
            shell("{path_picard}picard CreateSequenceDictionary R={input} O={output}  2>{log} 1>{log} ")