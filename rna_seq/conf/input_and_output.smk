#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: BioPipe
# Script : input_and_output.py
# Author : Peng Jia
# Date   : 2020.11.26
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
path_data = "/home/pengjia/ProjectSnake/rna_test"
path_log = "/home/pengjia/ProjectSnake/rna_test/logs"
path_data = path_data.rstrip("/") + "/"
path_log = path_log.rstrip("/") + "/"

reference_version = "GRCh38"
path_genome = "/home/DATA/REFGENOMEDB/human/GRCh38.d1.vd1/genome/GRCh38.d1.vd1.fa"
path_gtf = "/home/DATA/REFGENOMEDB/human/GRCh38.d1.vd1/annotation/GRCh38.96/gtf/Homo_sapiens.GRCh38.96.addChr.gtf"
path_gff = "/home/DATA/REFGENOMEDB/human/GRCh38.d1.vd1/annotation/GRCh38.96/gff/Homo_sapiens.GRCh38.96.addChr.gff"
path_STAR_index = path_data + "genome/star/index"

path_tophat2_index = path_data + "genome/tophat2/genome"
