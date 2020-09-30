# ======================================================================================================================
# Project: DNAseq
# Script : raw_read_qc_output.smk TODO check
# Author : Peng Jia
# Date   : 2020.07.21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
include: "raw_read_qc.smk"
path_raw_qc = []
path_qc_report = []
path_clean_reads = []
trim_method = config["pipe"]["trim"]

path_qc_report += expand(
    [path_data + "raw_read_qc/fastqc/{u.sample}/{u.sample}_{u.unit}_1_fastqc.zip"],
    u=caseinfo.itertuples()) + expand(
    [path_data + "raw_read_qc/fastqc/{u.sample}/{u.sample}_{u.unit}_2_fastqc.zip"],
    u=caseinfo.itertuples())
if config["pipe"]["trim"]!="passqc":
    path_qc_report +=expand(
    [path_data + "raw_read_qc/" + trim_method + "_fastqc2/{u.sample}/{u.sample}_{u.unit}_" +
     trim_method + "_fastqc2_1.zip"], u=caseinfo.itertuples()) + expand(
    [path_data + "raw_read_qc/" + trim_method + "_fastqc2/{u.sample}/{u.sample}_{u.unit}_" +
     trim_method + "_fastqc2_2.zip"], u=caseinfo.itertuples()) + expand(
    [path_data + "raw_read_qc/" + trim_method + "/{u.sample}/{u.sample}_{u.unit}_" +
     trim_method + ".json"],
    u=caseinfo.itertuples())
path_clean_reads = expand([path_data + "cleandata/" + trim_method +
                           "/{u.sample}/{u.sample}_{u.unit}_" +
                           trim_method + "_1.fq.gz"], u=caseinfo.itertuples()) + \
                   expand([path_data + "cleandata/" + trim_method +
                           "/{u.sample}/{u.sample}_{u.unit}_" +
                           trim_method + "_2.fq.gz"], u=caseinfo.itertuples())

path_raw_qc = path_qc_report + path_clean_reads
