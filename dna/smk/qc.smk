localrules: passqc

rule loadrawData:
    input:
         unpack(get_fastq)
    output:
          R1=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_1.fq.gz",
          R2=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_2.fq.gz"
    threads:
           config["threads"]["loadrawData"]
    log:
       path_log + "qc/loadRawdata/{case}/{sample}/{case}_{sample}_{unit}_loaddata.log"
    params:
          path=path_pigz
    wrapper:
           config["wrapper"] + "loadfq"

rule fastqc:
    input:
         R1=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_2.fq.gz"
    output:
          html1=path_data + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_1_fastqc.html",
          html2=path_data + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_2_fastqc.html",
          zip1=path_data + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_1_fastqc.zip",
          zip2=path_data + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_2_fastqc.zip"
    log:
       path_log + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_fastqc.log"
    threads: config["threads"]["fastqc"]
    benchmark:
             path_bm + "qc/fastqc/{case}/{sample}/{case}_{sample}_{unit}_fastqc.tsv"
    params:
          path=path_fastqc,
          extra=""
    wrapper:
           config["wrapper"] + "fastqc"

rule fastp:
    input:
         R1=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_2.fq.gz"
    output:
          R1=path_data + "cleandata/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp_1.fq.gz",
          R2=path_data + "cleandata/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp_2.fq.gz",
          html=path_data + "qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.html",
          json=path_data + "qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.json"
    log:
       path_log + "qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.log"
    threads: config["threads"]["fastp"]
    params:
          path=path_fastp,
          extra=""
    benchmark:
             path_bm + "qc/fastp/{case}/{sample}/{case}_{sample}_{unit}_fastp.tsv"
    wrapper:
           config["wrapper"] + "fastp"

rule passqc:
    input:
         R1=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{case}/{sample}/{case}_{sample}_{unit}_2.fq.gz"
    output:
          R1=path_data + "cleandata/passqc/{case}/{sample}/{case}_{sample}_{unit}_passqc_1.fq.gz",
          R2=path_data + "cleandata/passqc/{case}/{sample}/{case}_{sample}_{unit}_passqc_2.fq.gz"
    shell:
           """
           ln -sr {input.R1} {output.R1}
           ln -sr {input.R2} {output.R2}
           slepp 60
           touch -h {output.R1}
           touch -h {output.R2}           
           """

rule fastqc2:
    input:
         R1=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
         R2=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz",
    output:
          html1=path_data + "qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_1.html",
          html2=path_data + "qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_2.html",
          zip1=path_data + "qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_1.zip",
          zip2=path_data + "qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2_2.zip"
    log:
       path_log+"qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2.log"
    threads: config["threads"]["fastqc"]
    params:
          path=path_fastqc,
          extra=""
    benchmark:
             path_bm+"qc/{qcpipe}_fastqc2/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_fastqc2.tsv"
    wrapper:
           config["wrapper"] + "fastqc"

rule multiqc:
    input:
         expand([path_data + "qc/fastqc/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_1_fastqc.zip"],
                u=caseinfo.itertuples()),
         expand([path_data + "qc/fastqc/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_2_fastqc.zip"],
                u=caseinfo.itertuples()),
         expand([path_data + "qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_1.zip"],
                u=caseinfo.itertuples()),
         expand([path_data + "qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_2.zip"],
                u=caseinfo.itertuples()),
         expand([path_data + "qc/fastp/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp.json"],
                u=caseinfo.itertuples())
    output:
          report(path_data + "qc/multiqc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
       path_log+"qc/multiqc/multiqc.log"
    params:
          path=path_multiqc,
          extra=""
    threads: config["threads"]["multiqc"]
    wrapper:
           config["wrapper"] + "multiqc"
