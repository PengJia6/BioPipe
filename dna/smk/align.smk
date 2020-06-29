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

rule Bwa:
    input:
         R1=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
         R2=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz",
         ref=path_genome,
         sindex=path_genome+".fai",
         bindex=path_genome+".bwt"
    output:
          path_data + "align/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa_sorted.bam",
    log:
       path_log + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa.logs"
    benchmark:
             path_bm + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}.bwa.tsv"
    params:
          index=path_genome,
          extra=get_read_group,
          sort="samtools",
          sort_order="coordinate",
          sort_extra=" -m 2G ",
          pathsamtools=path_samtools,
          pathbwa=path_bwa
    threads: 8
    wrapper:
           config["wrapper"] + "bwa/mem"

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
    threads: 8
    params:
          path=path_samtools,
          extra=" ",
    wrapper:
           config["wrapper"] + "samtools/merge"

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
    threads: config["threads"]["picard"]["rmDup"]
    params:
          path=path_picard,
          java_opts=" -Xmx10G ",
          extra=""
    wrapper:
           config["wrapper"] + "picard/markduplicates"

rule NoneMarkDup:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_noneMarkDup.bam",
    shell:
         """
         ln -sr {input.bam} {output.bam}
         sleep 60
         touch -h {output.bam}
         """

rule MarkDupWithBiobambam:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.bam",
          metrics=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.Metrics",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bioMarkDup.tsv"
    threads: config["threads"]["biobambam"]["rmDup"]
    params:
          path=path_biobambam,
          extra=" tmpfile=../data/align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioMarkDup_tmp ",

    run:
        shell("{path_biobambam}bammarkduplicates I={input} O={output.bam} M={output.metrics} "
              "markthreads={threads} {params.extra} 1>{log} 2>{log}")
           # config["wrapper"] + "biobambam/bammarkduplicates"

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
          sleep 60
          touch -h {output.bam}
          touch -h {output.bai}
          """

rule GATKReAlignPre:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=path_genome,
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
    output:
          path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.intervals"
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_reAlign.tsv"
    threads: config["threads"]["gatk"]["preReAlign"]
    params:
          path=path_gatk3,
          extra="",  # optional
          java_opts=" -Xmx10G ",
    wrapper:
           config["wrapper"] + "/gatk3/realignertargetcreator"

rule GATKReAlign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=path_genome,
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
          path=path_gatk3,
          do=r"reAlign",  # reAlign or noAlign
          extra="",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["applyReAlign"]
    wrapper:
           config["wrapper"] + "/gatk3/indelrealigner"

rule BamIndex:
    input:
         "{prefix}.bam"
    output:
          "{prefix}.bam.bai"
    threads:config["threads"]["samtools"]
    params:
          path=path_samtools,
          extra=""
    wrapper:
           config["wrapper"] + "samtools/index"

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
         sleep 60
         touch -h {output.bam}
         touch -h {output.bai}
         """

rule GATKBQSRPre:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         ref=path_genome,
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
    output:
          target_intervals=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.tsv"
    params:
          path=path_gatk,
          extra=" --tmp-dir ../data/align/{case}/{sample}/ ",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["preBQSR"]
    wrapper:
           config["wrapper"] + "/gatk/baserecalibrator"

rule GATKBQSR:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         ref=path_genome,
         target_intervals=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals"
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.bam",
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.tsv"
    params:
          path=path_gatk,
          extra=" --tmp-dir ../data/align/{case}/{sample}/ ",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["applyBQSR"]
    wrapper:
           config["wrapper"] + "/gatk/applyBQSR"
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
         sleep 60
         touch -h {output.bam}
         touch -h {output.bai}
         """

rule GATKLeftAlign:
    input:
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}.bam.bai",
         ref=path_genome
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_leftAlign.bam"
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_leftAlign.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_leftAlign.tsv"
    params:
          path=path_gatk,
          extra=" --tmp-dir ../data/align/{case}/{sample}/ ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["leftAlign"]
    wrapper:
           config["wrapper"] + "/gatk/leftAlignIndels"
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
             sleep 60
             touch -h {output.bam}
             touch -h {output.bai}
             """

rule GATKFixMate:
    input:
         ref=path_genome,
         bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam",
         bai=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}.bam.bai"
    output:
          bam=path_data + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_fixMate.bam"
    log:
       path_log + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_fixMate.logs"
    benchmark:
             path_bm + "align/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_{bqsr}_{leftAlign}_fixMate.tsv"
    params:
          path=path_gatk,
          extra=" --TMP_DIR ../data/align/{case}/{sample}/ ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["fixMate"]
    wrapper:
           config["wrapper"] + "/gatk/fixMateInformation"

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
         sleep 60
         touch -h {output.bam}
         touch -h {output.bai}
         """
