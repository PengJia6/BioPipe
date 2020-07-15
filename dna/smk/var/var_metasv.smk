# ==============================================================================
# Project: DNAseq
# Script : var_metasv.smk
# Author : Peng Jia
# Date   : 2020.07.14
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Variants detection and filtering using metasv
# Homepage of metasv: https://github.com/bioinform/metasv
# [AGE] Abyzov,A. and Gerstein,M. (2011) AGE: defining breakpoints of genomic structural variants at single-nucleotide resolution, through optimal alignments with gap excision. Bioinformatics, 27, 595–603.
# [BreakDancer] Chen,K. et al. (2009) BreakDancer: an algorithm for high-resolution mapping of genomic structural variation. Nat. Methods, 6, 677–681.
# [BreakSeq2] Abyzov,A. et al. (2015) Analysis of deletion breakpoints from 1,092 humans reveals details of mutation mechanisms. Nat. Commun., 6, 7256.
# [BreakSeq] Lam,H.Y. et al. (2010) Nucleotide-resolution analysis of structural variants using BreakSeq and a breakpoint library. Nat. Biotechnol., 28, 47–55.
# [CNVnator] Abyzov,A. et al. (2011) CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing. Genome Res., 21, 974–984.
# [Pindel] Ye,K. et al. (2009) Pindel: a pattern growth approach to detect break points of large deletions and medium sized insertions from paired-end short reads. Bioinformatics, 25, 2865–2871.
# [SPAdes] Ye,K. et al. (2009) Bankevich,A. et al. (2012) SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J. Comput. Biol., 19, 455–477.
# ==============================================================================

# localrules: a
# ruleorder: a > b


# rules: MetaSV
# description: Variants Calling by MetaSV
# input: bam file
# output: raw vcf file
rule MetaSV:
    input:
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          workdir=directory(path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv_work"),
          outdir=directory(path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv_out"),
          cnvnator=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.cnvnator",
          breaddancer=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.breakdancer",
          breakseq=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.breakseq",
          pindel_D=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.pindel_D",
          pindel_LI=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.pindel_LI",
          pindel_SI=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.pindel_SI",
          pindel_TD=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.pindel_TD",
          pindel_INV=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}_metasv.pindel_INV",
          vcf=path_data + "germlineVar/metasv/perSample/{bam_sample}/{bam_sample}.metasv.raw.vcf.gz"
    log:
       path_log + "gremlineVar/metasv/perSample/{bam_sample}/{bam_sample}.metasv.logs"
    benchmark:
             path_bm + "gremlineVar/metasv/perSample/{bam_sample}/{bam_sample}.metasv.logs"

    threads: config["threads"]["MetaSV"]
    params:
          extra="",
          min_ins_support=2,
          max_ins_intervals=500000,
          isize_mean=500,
          isize_sd=150,

    run:
        shell("{path_metasv}run_metasv.py "
              "--reference {input.ref} "
              "--boost_sc "
              "--breakdancer_native {output.breaddancer} "
              "--breakseq_native {output.breakseq} "
              "--cnvnator_native {output.cnvnator} "
              "--pindel_native {output.pindel_D} {output.pindel_LI} {output.pindel_SI} {output.pindel_TD} {output.pindel_INV} "
              "--sample {wildcards.bam_sample} --bam {input.bam} "
              "--spades {path_spades}spades.py --age {path_age}/age_align "
              "--num_threads {threads} --workdir {output.workdir} "
              "--outdir {output.outdir} "
              "--min_support_ins {params.min_ins_support} "
              "--max_ins_intervals {params.max_ins_intervals} "
              "--isize_mean {params.isize_mean} "
              "--isize_sd {params.isize_sd} 2>{log} 1>{log}")

# # rules: TODO
# # description: TODO
# # input: TODO
# # output: TODO
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
#     shell:
#          """
#          command 1
#          command 2
#          """
