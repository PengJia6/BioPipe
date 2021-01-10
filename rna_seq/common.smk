# ======================================================================================================================
# Project: BioPipe
# Script : common.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
configfile: "conf/config.yaml"
include: "conf/software.smk"
include: "conf/input_and_output.smk"
include: "genome_database.smk"
include: "sample.smk"
include: "qc/loaddata.smk"
include: "qc/fastp.smk"
include: "align/star.smk"
include: "gene_fusion/STAR_fusion.smk"

wildcard_constraints:
                    qcpipe="fastp|trim|passqc",
                    aligner="STAR|bowtie",
                    markdup="bioMarkDup|picardMarkDup|noneMarkDup",
                    R="1|2",
                    sample="|".join(caseinfo["sample"]),
                    unit="|".join(caseinfo["unit"]),
                    LB="|".join(caseinfo["LB"]),

read_lens = list(set(list(map(int, (caseinfo["read_len"] - 1)))))
# print(read_lens)
output_files_genomes = []
path_align = []
path_fusion = []
sample_list = list(set([i[0] for i in caseinfo.index]))
sample_list=["1443727-N"]
if "STAR" in config["align"]:
    for read_len in read_lens:
        output_files_genomes.append(path_STAR_index + "_{read_len}".format(read_len=read_len))
    qc_pipe = config["qc"]
    for sample in sample_list:
        path_align.append(path_data + "align/STAR/{sample}/{sample}_{qcpipe}_{aligner}.bam".format(sample=sample,
                                                                                                   qcpipe=qc_pipe,
                                                                                                   aligner="STAR"))
# print(path_align)

if "STARfusion" in config["fusion"]:
    qc_pipe = config["qc"]
    for sample in sample_list:
        path_fusion.append(path_data +  "fusion/STARfusion/{sample}/{sample}_{qcpipe}_STAR.STARfusion".format(sample=sample,
                                                                                               qcpipe=qc_pipe))
if "tophat2" in config["align"]:
    path_align.append(path_tophat2_index)

