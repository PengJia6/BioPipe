localrules: passqc

import os

rule loadrawData:
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
             path_bm + "raw_read_qc/loadRawdata/{sample}/{sample}_{unit}_loaddata.tsv"
    params:
          ""
    run:
        inR1 = input.R1
        inR2 = input.R2
        if inR1.endswith("gz"):
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

rule fastqc:
    input:
         R1=path_data + "rawdata/{sample}/{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{sample}/{sample}_{unit}_2.fq.gz"
    output:
          tmp_dir=directory(path_data + "raw_read_qc/fastqc/{sample}/{sample}_{unit}/"),
          html1=path_data + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_1_fastqc.html",
          html2=path_data + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_2_fastqc.html",
          zip1=path_data + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_1_fastqc.zip",
          zip2=path_data + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_2_fastqc.zip"
    log:
       path_log + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_fastqc.log"
    threads: config["threads"]["fastqc"]
    benchmark:
             path_bm + "raw_read_qc/fastqc/{sample}/{sample}_{unit}_fastqc.tsv"
    params:
          path=path_fastqc,
          extra=""
    run:
        if not os.path.exists(output.tmp_dir):
            shell("mkdir {output.tmp_dir}")

        shell("{path_fastqc}fastqc {params.extra}  -t {threads} "
              "--outdir {output.tmp_dir} {input.R1} {input.R2}"
              " 2>>{log} 1>>{log}")


        def basename_without_ext(file_path):
            """Returns basename of file path, without the file extension."""
            base = os.path.basename(file_path)
            split_ind = 2 if base.endswith(".gz") else 1
            base = ".".join(base.split(".")[:-split_ind])
            return base
        output_base = basename_without_ext(input.R1)
        html_pathR1 = os.path.join(output.tmp_dir, output_base + "_fastqc.html")
        zip_pathR1 = os.path.join(output.tmp_dir, output_base + "_fastqc.zip")
        if output.html1 != html_pathR1:
            shell("mv {html_pathR1} {output.html1}")
        if output.zip1 != zip_pathR1:
            shell("mv {zip_pathR1} {output.zip1}")
        output_base = basename_without_ext(input.R2)
        zip_pathR2 = os.path.join(output.tmp_dir, output_base + "_fastqc.zip")
        html_pathR2 = os.path.join(output.tmp_dir, output_base + "_fastqc.html")
        if output.html2 != html_pathR2:
            shell("mv {html_pathR2} {output.html2}")
        if output.zip2 != zip_pathR2:
            shell("mv {zip_pathR2} {output.zip2}")

rule fastp:
    input:
         R1=path_data + "rawdata/{sample}/{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{sample}/{sample}_{unit}_2.fq.gz"
    output:
          R1=path_data + "cleandata/fastp/{sample}/{sample}_{unit}_fastp_1.fq.gz",
          R2=path_data + "cleandata/fastp/{sample}/{sample}_{unit}_fastp_2.fq.gz",
          html=path_data + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.html",
          json=path_data + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.json"
    log:
       path_log + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.log"
    threads: config["threads"]["fastp"]
    params:
          path=path_fastp,
          extra=""
    benchmark:
             path_bm + "raw_read_qc/fastp/{sample}/{sample}_{unit}_fastp.tsv"
    run:
        shell("{path_fastp}fastp --thread {threads} {params.extra} "
              "--in1 {input.R1} --in2 {input.R2} "
              "--out1 {output.R1} --out2 {output.R2} "
              "--html {output.html} --json {output.json} "
              "1>{log} 2>{log}")
# wrapper:
#        config["wrapper"] + "fastp"

rule passqc:
    input:
         R1=path_data + "rawdata/{sample}/{sample}_{unit}_1.fq.gz",
         R2=path_data + "rawdata/{sample}/{sample}_{unit}_2.fq.gz"
    output:
          R1=path_data + "cleandata/passqc/{sample}/{sample}_{unit}_passqc_1.fq.gz",
          R2=path_data + "cleandata/passqc/{sample}/{sample}_{unit}_passqc_2.fq.gz"
    shell:
         """
         ln -sr {input.R1} {output.R1}
         ln -sr {input.R2} {output.R2}
         sleep 1
         touch -h {output.R1}
         touch -h {output.R2}           
         """

rule fastqc2:
    input:
         R1=path_data + "cleandata/{qcpipe}/{sample}/{sample}_{unit}_{qcpipe}_1.fq.gz",
         R2=path_data + "cleandata/{qcpipe}/{sample}/{sample}_{unit}_{qcpipe}_2.fq.gz",
    output:
          tmp_dir=directory(path_data + "raw_read_qc/fastqc2/{sample}/{sample}_{unit}_{qcpipe}/"),
          html1=path_data + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2_1.html",
          html2=path_data + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2_2.html",
          zip1=path_data + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2_1.zip",
          zip2=path_data + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2_2.zip"
    log:
       path_log + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2.log"
    threads: config["threads"]["fastqc"]
    params:
          path=path_fastqc,
          extra=""
    benchmark:
             path_bm + "raw_read_qc/{qcpipe}_fastqc2/{sample}/{sample}_{unit}_{qcpipe}_fastqc2.tsv"
    run:
        if not os.path.exists(output.tmp_dir):
            shell("mkdir {output.tmp_dir}")

        shell("{path_fastqc}fastqc {params.extra}  -t {threads} "
              "--outdir {output.tmp_dir} {input.R1} {input.R2}"
              " 2>{log} 1>{log}")


        def basename_without_ext(file_path):
            """Returns basename of file path, without the file extension."""
            base = os.path.basename(file_path)
            split_ind = 2 if base.endswith(".gz") else 1
            base = ".".join(base.split(".")[:-split_ind])
            return base


        output_base = basename_without_ext(input.R1)
        html_pathR1 = os.path.join(output.tmp_dir, output_base + "_fastqc.html")
        zip_pathR1 = os.path.join(output.tmp_dir, output_base + "_fastqc.zip")
        if output.html1 != html_pathR1:
            shell("mv {html_pathR1} {output.html1}")
        if output.zip1 != zip_pathR1:
            shell("mv {zip_pathR1} {output.zip1}")
        output_base = basename_without_ext(input.R2)
        zip_pathR2 = os.path.join(output.tmp_dir, output_base + "_fastqc.zip")
        html_pathR2 = os.path.join(output.tmp_dir, output_base + "_fastqc.html")
        if output.html2 != html_pathR2:
            shell("mv {html_pathR2} {output.html2}")
        if output.zip2 != zip_pathR2:
            shell("mv {zip_pathR2} {output.zip2}")
#
# rule multiqc:
#     input:
#          expand([path_data + "raw_read_qc/fastqc/{u.sample}/{u.case}_{u.sample}_{u.unit}_1_fastqc.zip"],
#                 u=caseinfo.itertuples()),
#          expand([path_data + "raw_read_qc/fastqc/{u.sample}/{u.case}_{u.sample}_{u.unit}_2_fastqc.zip"],
#                 u=caseinfo.itertuples()),
#          expand([
#              path_data + "raw_read_qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_1.zip"],
#              u=caseinfo.itertuples()),
#          expand([
#              path_data + "raw_read_qc/fastp_fastqc2/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp_fastqc2_2.zip"],
#              u=caseinfo.itertuples()),
#          expand([path_data + "raw_read_qc/fastp/{u.case}/{u.sample}/{u.case}_{u.sample}_{u.unit}_fastp.json"],
#                 u=caseinfo.itertuples())
#     output:
#           report(path_data + "raw_read_qc/multiqc/multiqc.html", caption="../report/multiqc.rst",
#                  category="Quality control")
#     log:
#        path_log + "raw_read_qc/multiqc/multiqc.log"
#     params:
#           path=path_multiqc,
#           extra=""
#     threads: config["threads"]["multiqc"]
#     # wrapper:
#     #        config["wrapper"] + "multiqc"
#     run:
#         input_dirs = set(os.path.dirname(fp) for fp in input)
#         output_dir = os.path.dirname(output[0])
#         output_name = os.path.basename(output[0])
#         shell("{path_multiqc}multiqc"
#               " {params.extra}"
#               " --force  -o {output_dir}"
#               " -n {output_name}"
#               " {input_dirs}"
#               " 2>{log} 1>{log}")
