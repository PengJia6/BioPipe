rule Bwa:
    input:
         R1=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
         R2=path_data + "cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz"
    output:
          path_data + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa_sorted.bam",
    log:
       path_log+"align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa.logs"
    benchmark:
             path_bm + "align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}.bwa.tsv"
    params:
          index=config["ref"]["genome"],
          extra=get_read_group,
          sort="samtools",
          sort_order="coordinate",
          sort_extra=" -m 2G ",
          pathsamtools=config["mainEnv"],
          pathbwa=config["mainEnv"]
    threads: 8
    wrapper:
           config["wrapper"] + "bwa/mem"

rule MergeBam:
    input:
         get_sample_bam
         #        expand("{path_prefix}align/{{aligner}}/{{case}}/{{sample}}/{{case}}_{{sample}}_{unit}_{{qcpipe}}_{{aligner}}_sorted.bam", unit=caseinfo.loc[({case},{sample}),"unit"],path_prefix=[path_data])
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.tsv"
    threads: 8
    params:
          path=config["mainEnv"],
          extra=" ",
    wrapper:
           config["wrapper"] + "samtools/merge"

rule MarkDupWithPicard:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_pcdRmDup.bam",
          metrics=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_pcdRmDup.metrics",
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardRmdup.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardRmdup.tsv"
    threads: config["threads"]["picard"]["rmDup"]
    params:
          path=config["mainEnv"],
          java_opts=" -Xmx10G ",
          extra=""
    wrapper:
           config["wrapper"] + "gatk/markduplicates"

rule MarkDupWithBiobambam:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup.bam",
          metrics=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDupMetrics",
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.BioRmdup.tsv"
    threads: config["threads"]["biobambam"]["rmDup"]
    params:
          path=config["mainEnv"],
          extra=" tmpfile=../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup_tmp ",
    wrapper:
           config["wrapper"] + "biobambam/bammarkduplicates"

rule GATKReAlignPre:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=config["ref"]["genome"],
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
    output:
          path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.intervals"
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.tsv"
    threads: config["threads"]["gatk"]["preReAlign"]
    params:
          path=config["mainEnv"],
          do=r"{realign}",  # realign or no align
          extra="",  # optional
          java_opts=" -Xmx10G ",
    wrapper:
           config["wrapper"] + "/gatk3/realignertargetcreator"

rule GATKReAlign:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}.bam.bai",
         ref=config["ref"]["genome"],
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
         target_intervals=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.intervals"
    output:
          path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam"
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}-apply.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}-apply.tsv"
    params:
          path=config["mainEnv"],
          do=r"{realign}",  # reAlign or noAlign
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
          path=config["mainEnv"],
          extra=""
    wrapper:
           config["wrapper"] + "samtools/index"

rule GATKBQSRPre:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         ref=config["ref"]["genome"],
         known=[config["ref"]["1KGp1snp"], config["ref"]["dbsnp"], config["ref"]["1KGomni"], config["ref"]["hapmap"],
                config["ref"]["mills1KG"]],
    output:
          target_intervals=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals",
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_preBQSR.tsv"
    params:
          path=config["mainEnv"],
          extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ ",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["preBQSR"]
    wrapper:
           config["wrapper"] + "/gatk/baserecalibrator"

rule GATKBQSR:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}.bam",
         ref=config["ref"]["genome"],
         target_intervals=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.intervals"
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.bam",
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR.tsv"
    params:
          path=config["mainEnv"],
          extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ ",  # optional
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["applyBQSR"]
    wrapper:
           config["wrapper"] + "/gatk/applyBQSR"
rule GATKLeftAlign:
    input:
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR.bam.bai",
         ref=config["ref"]["genome"]
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR_leftAlign.bam"
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR_leftAlign.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR_leftAlign.tsv"
    params:
          path=config["mainEnv"],
          extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["leftAlign"]
    wrapper:
           config["wrapper"] + "/gatk/leftAlignIndels"
rule GATKFixMate:
    input:
         ref=config["ref"]["genome"],
         bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR_leftAlign.bam",
         bai=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR_leftAlign.bam.bai"
    output:
          bam=path_data + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_BQSR_leftAlign_fixMate.bam"
    log:
       path_log+"align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR_leftAlign_fixMate.logs"
    benchmark:
             path_bm + "align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{markdup}_{realign}_applyBQSR_leftAlign_fixMate.tsv"
    params:
          path=config["mainEnv"],
          extra=" --TMP_DIR ../data/align/{aligner}/{case}/{sample}/ ",
          java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["fixMate"]
    wrapper:
           config["wrapper"] + "/gatk/fixMateInformation"
