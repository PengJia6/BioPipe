# ======================================================================================================================
# Project: BioPipe
# Script : software.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
default = "/home/pengjia/miniconda3/envs/RNAseq/bin/"

default = default.rstrip("/") + "/"
path_star = default + "STAR "
path_star_fusion = default + "STAR-Fusion "  # offline installation
path_fastp = default + "fastp "

path_tophat2 = "/home/pengjia/miniconda3/envs/tophat/bin/tophat2"
