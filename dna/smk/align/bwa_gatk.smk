localrules: LoadHQbam
localrules: NoneBQSR
localrules: NoneLeftALign
localrules: NoneFixMate
localrules: NoneRealign

ruleorder: LoadHQbam > BamIndex
ruleorder: NoneBQSR > BamIndex
ruleorder: NoneLeftALign > BamIndex
ruleorder: NoneFixMate > BamIndex
ruleorder: NoneRealign > BamIndex

import os

# rules: Bwa
# description: Bwa and sort
# input: clean.fq.gz
# output: sorted bam
# check: PASS
rule Bwa:
    input:
         R1=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
         R2=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz",
         ref=path_genome,
         sindex=path_genome + ".fai",
         bindex=path_genome + ".bwt"
    output:
          path_data + "align/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa_sorted.bam",
    log:
       path_log + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa.logs"
    benchmark:
             path_bm + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}.bwa.tsv"
    params:
          extra=get_read_group,
          sort_extra=" -m 2G ",
          bwa_extra=""
    threads: config["threads"]["Bwa"]
    run:
        shell("{path_bwa}bwa mem -M {params.extra} -t {threads} {input.ref} {input.R1} {input.R2} | "
              "{path_samtools}samtools view -Shb -@ {threads} | "
              "{path_samtools}samtools sort -@ {threads} {params.sort_extra} -T {output}_tmp -o {output} -O BAM "
              "1>{log} 2>{log} ")

# rules: MergeBam
# description: merge bam for different Lane and lib
# input: bam
# output:  merged bam
# check: PASS
rule MergeBam:
    input:
         get_sample_bam
         #        expand("{path_prefix}align/{{aligner}}/{{case}}/{{sample}}/{{case}}_{{sample}}_{unit}_{{qcpipe}}_{{aligner}}_sorted.bam", unit=caseinfo.loc[({case},{sample}),"unit"],path_prefix=[path_data])
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.tsv"
    threads: config["threads"]["MergeBam"]
    params:
          # path=path_samtools,
          extra=" ",
    run:
        if len(input) < 2:
            shell("ln -sr {input} {output.bam} ")
            shell("sleep 1")
            shell("touch -h {output.bam} 1>>{log} 2>>{log}")
            shell("echo only one bam, make soft link 1>>{log} 2>>{log}")
        else:
            shell("{path_samtools}samtools merge -@ {threads} {params.extra} "
                  "{output.bam} {input} 2>{log} 1>{log}")

# wrapper:
#        config["wrapper"] + "samtools/merge"
# rules: MarkDupWithPicard
# description: mark duplicated reads with picard
# input: bam
# output:  marked bam
# check: TODO
rule MarkDupWithPicard:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardMarkDup.bam",
          metrics=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardMarkDup.metrics",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardMarkDup.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardMarkDup.tsv"
    threads: config["threads"]["MarkDupWithPicard"]
    params:
          java_opts=" -Xmx10G ",
          extra=""
    run:
        shell("{path_picard}picard  MarkDuplicates {params.extra} "
              "I={input} "
              "O={output.bam} M={output.metrics} "
              " 1>{log} 2>{log}")

# wrapper:
#        config["wrapper"] + "picard/markduplicates"

# rules: NoneMarkDup
# description: do not mark duplicated reads
# input: bam
# output:  bam
# check: PASS
rule NoneMarkDup:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_noneMarkDup.bam",
    shell:
         """
         ln -sr {input.bam} {output.bam}
         sleep 1
         touch -h {output.bam}
         """
# rules: MarkDupWithBiobambam
# description: mark duplicated reads with Biobambam
# input: bam
# output:  marked bam
# check: PASS
rule MarkDupWithBiobambam:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          tmp=directory(path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_biobam_tmp"),
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.bam",
          metrics=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.Metrics",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bioMarkDup.tsv"
    threads: config["threads"]["MarkDupWithBiobambam"]
    params:
          extra="",
    run:
        if not os.path.exists(output.tmp):
            shell("mkdir {output.tmp}")
        shell("{path_biobambam}bammarkduplicates2 I={input} O={output.bam} M={output.metrics} "
              "markthreads={threads}  tmpfile={output.tmp} 1>{log} 2>{log}")

# rules: NoneRealign
# description: do not realign reads loated in indels regions
# input: bam
# output:  bam
# check: PASS
rule NoneRealign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_noneReAlign.bam",
          bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.noneReAlign.bam.bai",
    shell:
         """
          ln -sr {input.bam} {output.bam}
          ln -sr {input.bai} {output.bai}
          sleep 1
          touch -h {output.bam}
          touch -h {output.bai}
          """

# rules: GATKReAlignPre
# description: realign reads located in indels regions with GATK
# input: bam
# output:  realigned bam
# check: TODO
rule GATKReAlignPre:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=path_genome,
         dict=path_dict,
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
    output:
          path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.intervals"
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.tsv"
    threads: config["threads"]["GATKReAlignPre"]
    params:
          extra="",  # optional
          java_opts=" -Xmx10G ",
          bed=""  # -L xxx.bed
    run:
        input_known_string = ""
        for known in input.known:
            input_known_string = input_known_string + " --known {}".format(known)
        shell("{path_gatk3}gatk3 {params.java_opts} -T RealignerTargetCreator"
              " -nt {threads} {params.extra} -I {input.bam} -R {input.ref}"
              " {input_known_string} {params.bed} -o {output} 2>{log} 1>{log}")

# wrapper:
#        config["wrapper"] + "/gatk3/realignertargetcreator"

rule GATKReAlign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=path_genome,
         dict=path_dict,
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
         target_intervals=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.intervals"
    output:
          path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.bam"
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign-apply.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign-apply.tsv"
    params:
          extra="",  # optional
          bed="",  # bed=""xx.bed
          java_opts=" -Xmx10G ",
    threads: config["threads"]["GATKReAlign"]
    run:
        input_known_string = ""
        for known in input.known:
            input_known_string = input_known_string + " -known {}".format(known)
        shell("{path_gatk3}gatk3 {params.java_opts} -T IndelRealigner"
              " {params.extra} -I {input.bam}  -R {input.ref}"
              " {input_known_string}  {params.bed} -o {output}"
              " --targetIntervals {input.target_intervals}"
              " 1>{log} 2>{log}")
# wrapper:
#        config["wrapper"] + "/gatk3/indelrealigner"


# rules: BamIndex
# description: bam index
# input: bam
# output:  bam.bai
# check: PASS
rule BamIndex:
    input:
         "{prefix}.bam"
    output:
          "{prefix}.bam.bai"
    threads:config["threads"]["BamIndex"]
    params:
          path=path_samtools,
          extra=""
    run:
        if os.path.exists(wildcards.prefix + ".bai"):
            shell("ln -sr {wildcards.prefix}.bai {output}")
        else:
            shell("{path_samtools}samtools index -@ {threads} {input} ")
# rules: NoneBQSR
# description: do not using BQSR model
# input: bam
# output:  bam
# check: PASS
rule NoneBQSR:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam.bai"
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_noneBQSR.bam",
          bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_noneBQSR.bam.bai"
    shell:
         """
         ln -sr {input.bam} {output.bam}
         ln -sr {input.bai} {output.bai}
         sleep 1
         touch -h {output.bam}
         touch -h {output.bai}
         """

# rules: GATKBQSRPre
# description: BQSR
# input: bam
# output:  bam
# check: TODO
rule GATKBQSRPre:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam.bai",
         ref=path_genome,
         dict=path_dict,
         # known=[config["ref"]["1KGomni"], config["ref"]["hapmap"],
         #        ],
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"]],
         # known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
         #             config["ref"]["mills1KG"]],
    output:
          tmp=directory(
              path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSRPre_tmp"),
          target_intervals=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.tsv"
    params:
          # path=path_gatk,
          extra="",
          # extra=" --tmp-dir ../data/align/{case}/{sample}/ ",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["GATKBQSRPre"]
    run:
        input_known_string = ""
        for known in input.known:
            input_known_string = input_known_string + " --known-sites {} ".format(known)
        if not os.path.exists(output.tmp):
            shell("mkdir {output.tmp}")
        shell("{path_gatk}gatk --java-options {params.java_opts} BaseRecalibrator {params.extra} "
              "--tmp-dir {output.tmp} -R {input.ref} -I {input.bam} "
              "-O {output.target_intervals} {input_known_string} "
              "1>{log} 2>{log}")
# wrapper:
#        config["wrapper"] + "/gatk/baserecalibrator"

rule GATKBQSR:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam.bai",
         ref=path_genome,
         target_intervals=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals"
    output:
          tmp=directory(
              path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR_tmp"),
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.bam",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.tsv"
    params:
          extra="",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["GATKBQSR"]
    run:
        if not os.path.exists(output.tmp):
            shell("mkdir {output.tmp}")
        shell("{path_gatk}gatk --java-options {params.java_opts} ApplyBQSR -R {input.ref} -I {input.bam} "
              "--bqsr-recal-file {input.target_intervals}  "
              "-O {output.bam} --tmp-dir {output.tmp} 1>{log} 2>{log}")

# rules: NoneLeftALign
# description: do not using left realign model
# input: bam
# output:  bam
# check: PASS
rule NoneLeftALign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam.bai",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_noneLeftAlign.bam",
          bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_noneLeftAlign.bam.bai"
    shell:
         """
         ln -sr {input.bam} {output.bam}
         ln -sr {input.bai} {output.bai}
         sleep 1
         touch -h {output.bam}
         touch -h {output.bai}
         """

# rules: GATKLeftAlign
# description: GATK left Align
# input: bam
# output:  bam
# check: TODO
rule GATKLeftAlign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam.bai",
         ref=path_genome
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_LeftAlign.bam",
          tmp=directory(
              path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_LeftAlign_tmp"),
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_leftAlign.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_leftAlign.tsv"
    params:
          extra=" ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["GATKLeftAlign"]
    run:
        if not os.path.exists(output.tmp):
            shell("mkdir {output.tmp}")
        shell("{path_gatk3}gatk3 {params.java_opts} -T LeftAlignIndels "
              " -R {input.ref} -I {input.bam}  "
              " -o {output.bam} 1>{log} 2>{log}")
# shell("{path_gatk}gatk --java-options   {params.java_opts} LeftAlignIndels "
#       " -R {input.ref} -I {input.bam}  "
#       " -O {output.bam} 1>{log} 2>{log}")


# rules: NoneFixMate
# description: do not fixmate
# input: bam
# output:  bam
# check: PASS
rule NoneFixMate:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam.bai"
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_noneFixMate.bam",
          bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_noneFixMate.bam.bai"
    shell:
         """
             ln -sr {input.bam} {output.bam}
             ln -sr {input.bai} {output.bai}
             sleep 1
             touch -h {output.bam}
             touch -h {output.bai}
             """
# rules: GATKFixMate
# description: Fix mate using gatk
# input: bam
# output:  bam
# check: TODO
rule GATKFixMate:
    input:
         ref=path_genome,
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam.bai"
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_FixMate.bam",
          tmp=directory(
              path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_FixMate_tmp")
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_fixMate.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_fixMate.tsv"
    params:
          extra="",
          # extra=" --TMP_DIR ../data/align/{case}/{sample}/ ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["GATKFixMate"]
    run:
        if not os.path.exists(output.tmp):
            shell("mkdir {output.tmp}")

        shell("{path_gatk}gatk --java-options {params.java_opts} FixMateInformation {params.extra} "
              " -R {input.ref} -I {input.bam} -SO coordinate "
              " -O {output.bam} 2>{log} 1>{log}")
# wrapper:
#        config["wrapper"] + "/gatk/fixMateInformation"


rule LoadHQbam:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_" + bam_suffix + ".bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_" + bam_suffix + ".bam.bai"
    output:
          bam=path_data + "HQbam/{case}_{sample}.bam",
          bai=path_data + "HQbam/{case}_{sample}.bam.bai",
          logs=path_data + "HQbam/changeLogs/{case}_{sample}.logs"

    shell:
         """
         echo `readlink -f {input.bam}` " ->" {output.bam} > {output.logs}
         ln -sr {input.bam} {output.bam}
         ln -sr {input.bai} {output.bai}
         sleep 1
         touch -h {output.bam}
         touch -h {output.bai}
         """
