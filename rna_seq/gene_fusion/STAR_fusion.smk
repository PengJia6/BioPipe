# ======================================================================================================================
# Project: BioPipe
# Script : STAR_fusion.smk
# Author : Peng Jia
# Date   : 2020.11.29
# Email  : pengjia@stu.xjtu.edu.cn
# Description:
# ======================================================================================================================

# localrules: a
# ruleorder: a > b

# ======================================================================================================================
# rules: TODO
# description: TODO
# input: TODO
# output: TODO
star_fusion_lib = "/home/DATA/REFGENOMEDB/human/GRCh38.d1.vd1/STAR/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir"
rule STARfusion:
    input:
         juction=path_data + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}_Chimeric.out.junction"
    output:
          tag=path_data + "fusion/STARfusion/{sample}/{sample}_{qcpipe}_{aligner}.STARfusion",
    log:
       path_log + "fusion/{sample}/{sample}_{qcpipe}_{aligner}.STARfusion.logs"
    benchmark:
             path_log + "fusion/{sample}/{sample}_{qcpipe}_{aligner}.STARfusion.logs"
    threads: 8
    params:
          extra="",
    run:
        star_path = "/".join(path_star.rstrip("\n").split("/")[:-1])
        shell(
            "export PATH={star_path}:$PATH && {path_star_fusion} --genome_lib_dir {star_fusion_lib} -J {input.juction} "
            "--output_dir {output.tag}_dir 2>{log} 1>{log} ")
        shell("touch {output.tag}")
