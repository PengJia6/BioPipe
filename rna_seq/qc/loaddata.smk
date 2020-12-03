# ======================================================================================================================
# Project: BioPipe
# Script : loaddata.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.27
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

# localrules: a
# ruleorder: a > b
# ======================================================================================================================

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = caseinfo.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"R1": fastqs.fq1, "R2": fastqs.fq2}
    return {"R1": fastqs.fq1, "R2": fastqs.fq2}


def get_read_len(wildcards):
    """Get fastq files of given sample-unit."""
    read_len = caseinfo.loc[(wildcards.sample, wildcards.unit), "read_len"]
    return int(read_len)


def get_sample_bams(wildcards):
    units = list(set(caseinfo.loc[(wildcards.sample), "unit"]))
    aligner = wildcards.aligner
    # case = wildcards.case
    sample = wildcards.sample
    qcpipe = wildcards.qcpipe
    # path_data + "align/"+aligner+"/{sample}/{sample}_{unit}_{qc_pipe}_STAR
    # bam = path_data + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR_Aligned.sortedByCoord.out.bam",

    return [path_data + "align/" + aligner + "/" + sample + "/" + sample + "_" + str(
        i) + "_" + qcpipe + "_" + aligner + "_Aligned.sortedByCoord.out.bam" for i in units]

def get_sample_jucntion(wildcards):
    units = list(set(caseinfo.loc[(wildcards.sample), "unit"]))
    aligner = wildcards.aligner
    # case = wildcards.case
    sample = wildcards.sample
    qcpipe = wildcards.qcpipe
    # path_data + "align/"+aligner+"/{sample}/{sample}_{unit}_{qc_pipe}_STAR
    # bam = path_data + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR_Aligned.sortedByCoord.out.bam",
    print( [path_data + "align/" + aligner + "/" + sample + "/" + sample + "_" + str(
        i) + "_" + qcpipe + "_" + aligner + "_Chimeric.out.junction" for i in units])
    return [path_data + "align/" + aligner + "/" + sample + "/" + sample + "_" + str(
        i) + "_" + qcpipe + "_" + aligner + "_Chimeric.out.junction" for i in units]
rule loadRawData:
    input:
         unpack(get_fastq)
    output:
          R1=path_data + "rawdata/{sample}/{sample}_{unit}_1.fq.gz",
          R2=path_data + "rawdata/{sample}/{sample}_{unit}_2.fq.gz"
    threads:
           config["threads"]["loadrawData"]
    log:
       path_log + "raw_read_qc/loadRawdata/{sample}/{sample}_{unit}_loaddata.log"
    benchmark:
             path_log + "raw_read_qc/loadRawdata/{sample}/{sample}_{unit}_loaddata.tsv"
    params:
          ""
    run:
        inR1 = input.R1
        inR2 = input.R2
        if str(inR1).endswith("gz"):
            shell("ln -sr {input.R1} {output.R1} 2>>{log} 1>>{log}")
            shell("echo file is compressed, make soft link 1>>{log} 2>>{log}")
        else:
            shell("echo Compress file using pigz 1>>{log} 2>>{log}")

            shell("{path_pigz}pigz -p {threads} -c {input.R1} >{output.R1}")
        if inR2.endswith("gz"):
            shell("ln -sr {input.R2} {output.R2} 2>>{log} 1>>{log}")
            shell("echo file is compressed, make soft link 1>>{log} 2>>{log}")
        else:
            shell("echo Compress file using pigz 1>>{log} 2>>{log}")
            shell("{path_pigz}pigz -p {threads} -c {input.R2} >{output.R2}")
