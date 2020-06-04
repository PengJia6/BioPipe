rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html1="../data/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_1_fastqc.html",
        html2="../data/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_2_fastqc.html",
        zip1="../data/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_1_fastqc.zip",
        zip2="../data/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_2_fastqc.zip"
    log:
        "../logs/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_fastqc.log"
    threads: config["threads"]["fastqc"]
    benchmark: 
        "../benchmark/qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_fastqc.tsv"
    params: 
        path=config["qcEnv"],
        extra=""
    wrapper:
        config["wrapper"]+"fastqc"

rule fastp: 
    input: 
        unpack(get_fastq)
    output:
        R1="../data/cleandata/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp_1.fq.gz",
        R2="../data/cleandata/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp_2.fq.gz",
        html="../data/qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.html",
        json="../data/qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.json"
    log: 
        "../logs/qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.log" 
    threads: config["threads"]["fastp"] 
    params: 
        path=config["qcEnv"],
        extra=""
    benchmark:
        "../benchmark/qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.tsv" 
    wrapper:
        config["wrapper"]+"fastp"

rule passqc: 
    input: 
        unpack(get_fastq)
    output:
        R1="../data/cleandata/passqc/{case}/{sample}/{case}_{sample}_{unit}_passqc_1.fq.gz",
        R2="../data/cleandata/passqc/{case}/{sample}/{case}_{sample}_{unit}_passqc_2.fq.gz"
    threads: 1
    wrapper:
        config["wrapper"]+"passqc"        


rule fastqc2:
    input: 
        R1="../data/cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
        R2="../data/cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz",
    output:
        html1="../data/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_1.html",
        html2="../data/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_2.html",
        zip1="../data/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_1.zip",
        zip2="../data/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_2.zip"
    log:
        "../logs/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2.log"
    threads: config["threads"]["fastqc"]
    params: 
        path=config["qcEnv"],
        extra=""
    benchmark: 
        "../benchmark/qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2.tsv"
    wrapper:
        config["wrapper"]+"fastqc"


rule multiqc:
    input: 
        expand(["../data/qc/fastqc/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_1_fastqc.zip"], u=caseinfo.itertuples()),
        expand(["../data/qc/fastqc/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_2_fastqc.zip"], u=caseinfo.itertuples()),
        expand(["../data/qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_1.zip"], u=caseinfo.itertuples()),
        expand(["../data/qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_2.zip"], u=caseinfo.itertuples()),
        expand(["../data/qc/fastp/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp.json"], u=caseinfo.itertuples())
    output:
        report("../data/qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "../logs/qc/multiqc/multiqc.log"
    params:
        path=config["qcEnv"],
        extra=""
    threads: config["threads"]["multiqc"]
    wrapper:
        config["wrapper"]+"multiqc"
