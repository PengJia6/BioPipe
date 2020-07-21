# ======================================================================================================================
# Project: DNAseq
# Script : align_output.smk TODO check 
# Author : Peng Jia
# Date   : 2020.07.21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
include: "bwa_gatk.smk"

path_HQbamList = expand([path_data + "HQbam/{u.case}_{u.sample}.bam.bai"], u=caseinfo.itertuples())

