# ======================================================================================================================
# Project: DNAseq
# Script : var_cnvcaller.smk
# Author : Peng Jia
# Date   : 2020.08.01
# Email  : pengjia@stu.xjtu.edu.cn
# Description: CNVcaller
# ======================================================================================================================

# localrules: a
# ruleorder: a > b


CNVcaller_bin_size = 500
sex_chrom = "chrX"
CNVcaller_ref_index = "referenceDB." + str(CNVcaller_bin_size)
path_CNVcaller_ref_index = path_genome_path + CNVcaller_ref_index
path_orginal_CNVcaller_index = path_orginal_genome_path + CNVcaller_ref_index
ref_name = config["ref"]["genome"].split("/")[-1]
path_orginal_CNVcaller_kmer = path_orginal_genome_prefix + "cnvcaller.kmer.fa"
path_CNVcaller_dup = path_genome_prefix + "cnvcaller.link"

rule CNVcaler_IndexReference:
    input:
         path_genome,
    output:
          path_CNVcaller_ref_index,
    log:
       os.path.abspath(path_log + "genome/CNVcaler/reference/IndexReference.logs")
    benchmark:
             os.path.abspath(path_bm + "genome/CNVcaler/reference/IndexReference.bm")
    threads: config["threads"]["CNVcaler_IndexReference"]
    params:
          extra="",
    run:
        # inpath=os.path.abspath(input)
        if os.path.exists(path_orginal_CNVcaller_index):
            shell("ln -sr {path_orginal_CNVcaller_index} {output}")
            shell("sleep 5")
            shell("touch -h {output}")
        else:
            shell("""
                  cd {path_genome_path} 
                  {path_perl}perl {path_cnvcaller}bin/CNVReferenceDB.pl -w {CNVcaller_bin_size} {input} 2>{log} 1>{log}
                  """)

rule CNVcaller_make_kmer:
    input:
         path_genome
    output:
          path_genome_prefix + "cnvcaller.kmer.fa"
    log:
       path_log + "genome/CNVcaler/reference/IndexReference.logs"
    benchmark:
             os.path.abspath(path_bm + "genome/CNVcaler/reference/IndexReference.bm")
    threads: config["threads"]["CNVcaller_make_kmer"]
    run:
        if os.path.exists(path_orginal_genome_prefix + "cnvcaller.kmer.fa"):
            shell("ln -sr {path_orginal_genome_prefix}cnvcaller.kmer.fa {output}")
            shell("sleep 5")
            shell("touch -h {output}")
        else:
            shell("{path_python3}python {path_cnvcaller}bin/0.1.Kmer_Generate.py {input} "
                  "{CNVcaller_bin_size} {output} 2>{log} 1>{log}")

rule CNVcaller_make_dup_align:
    input:
         kmer=rules.CNVcaller_make_kmer.output,
         ref=path_genome,
    output:
          path_genome_prefix + "cnvcaller.kmer.aln"
    log:
       path_log + "genome/CNVcaler/reference/make_dup_align.logs"
    benchmark:
             path_bm + "genome/CNVcaler/reference/make_dup_aln.bm"
    threads: config["threads"]["CNVcaller_make_dup_align"]
    run:
        if os.path.exists(path_orginal_genome_prefix + "cnvcaller.kmer.aln"):
            shell("ln -sr {path_orginal_genome_prefix}cnvcaller.kmer.aln {output}")
            shell("sleep 5")
            shell("touch -h {output}")
        else:
            shell("{path_blasr}blasr {input.kmer} {input.ref} --sa {input.ref}.sa --out {output} "
                  " -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 --advanceHalf "
                  "--advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut "
                  "--bestn 10  2>>{log} 1>>{log}")

rule CNVcaller_make_dup:
    input:
         kmer=rules.CNVcaller_make_kmer.output,
         ref=path_genome,
         aln=rules.CNVcaller_make_dup_align.output
    output:
          path_genome_prefix + "cnvcaller.link"
    log:
       path_log + "genome/CNVcaler/reference/CNVcaller_make_dup.logs"
    benchmark:
             path_bm + "genome/CNVcaler/reference/CNVcaller_make_dup.bm"
    threads: config["threads"]["CNVcaller_make_dup"]
    run:
        if os.path.exists(path_orginal_genome_prefix + "cnvcaller.link"):
            shell("ln -sr {path_orginal_genome_prefix}cnvcaller.link {output}")
            shell("sleep 5")
            shell("touch -h {output}")
        else:
            # shell
            shell("")
            shell("{path_python3}python {path_cnvcaller}bin/0.2.Kmer_Link.py {input.aln} "
                  "{CNVcaller_bin_size} window.link 2>>{log} 1>>{log}")

# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO


# shell(
#     """
#     cd {path_CNVcaller_ref}
#    {path_perl} {path_cnvcaller}
#  """
# )
# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
rule Individual_RD_Processing:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
         dup=path_CNVcaller_ref_index,
         cnvcaller_index=path_CNVcaller_dup
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
