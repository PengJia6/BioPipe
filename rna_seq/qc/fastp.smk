# ======================================================================================================================
# Project: BioPipe
# Script : fastp.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.27
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
rule fastp:
    input:
         R1=path_data + "rawdata/{sample}/{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{sample}/{sample}_{unit}_2.fq.gz"
    output:
          R1=path_data + "cleandata/fastp/{sample}/{sample}_{unit}_fastp_1.fq.gz",
          R2=path_data + "cleandata/fastp/{sample}/{sample}_{unit}_fastp_2.fq.gz",
          html=path_data + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.html",
          json=path_data + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.json"

    threads: config["threads"]["fastp"]
    params:
          path=path_fastp,
          extra=""
    log:
       path_log + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.log"
    benchmark:
             path_log + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.tsv"
    run:
        shell("{path_fastp} --thread {threads} {params.extra} "
              "--in1 {input.R1} --in2 {input.R2} "
              "--out1 {output.R1} --out2 {output.R2} "
              "--html {output.html} --json {output.json} "
              "1>{log} 2>{log}")
