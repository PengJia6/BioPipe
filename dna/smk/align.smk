rule bwa: 
    input:
        R1="../data/cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_1.fq.gz",
        R2="../data/cleandata/{qcpipe}/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_2.fq.gz"
    output:
        "../data/align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa_sorted.bam",
    log:
        "../logs/align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}_bwa.logs"
    benchmark:
        "../benchmark/align/bwa/{case}/{sample}/{case}_{sample}_{unit}_{qcpipe}.bwa.tsv"
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
        config["wrapper"]+"bwa/mem"

rule bammerge: 
    input: 
        get_sample_bam
#        expand("../data/align/{{aligner}}/{{case}}/{{sample}}/{{case}}_{{sample}}_{unit}_{{qcpipe}}_{{aligner}}_sorted.bam", unit=caseinfo.loc[({case},{sample}),"unit"])
    output: 
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.merge.tsv"
    threads: 8
    params: 
        path=config["mainEnv"],
        extra=" ",
    wrapper:
        config["wrapper"]+"samtools/merge"


rule picardRmDup:
    input: 
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_pcdRmDup.bam",
        metrics= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_pcdRmDup.metrics",
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardRmdup.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_picardRmdup.tsv"
    threads: config["threads"]["picard"]["rmDup"]
    params:  
        path=config["mainEnv"],
        java_opts=" -Xmx10G ",
        extra=""
    wrapper:
        config["wrapper"]+"gatk/markduplicates"


rule bioRmDup:
    input: 
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.bam",
    output:
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup.bam",
        metrics= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDupMetrics",
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}.BioRmdup.tsv"
    threads: config["threads"]["biobambam"]["rmDup"]
    params: 
        path=config["mainEnv"],
        extra=" tmpfile=../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_bioRmDup_tmp ",
    wrapper:  
        config["wrapper"]+"biobambam/bammarkduplicates"

rule preReAlign:
    input:
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}.bam",
        bai= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}.bam.bai",
        ref=config["ref"]["genome"],
        known=[config["ref"]["1KGp1snp"],config["ref"]["dbsnp"],config["ref"]["1KGomni"],config["ref"]["hapmap"],config["ref"]["mills1KG"]],
    output:
        "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.intervals"
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.tsv"
    threads: config["threads"]["gatk"]["preReAlign"]
    params:
        path=config["mainEnv"],
        do=r"{realign}",   # realign or no align
        extra="",  # optional
        java_opts=" -Xmx10G ",
    wrapper:
           config["wrapper"]+"/gatk3/realignertargetcreator"

rule applyReAlign:
    input:
        bam= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}.bam",
        bai= "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}.bam.bai",
        ref=config["ref"]["genome"],
        known=[config["ref"]["1KGp1snp"],config["ref"]["dbsnp"],config["ref"]["1KGomni"],config["ref"]["hapmap"],config["ref"]["mills1KG"]],
        target_intervals="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.intervals"
    output:
        "../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.bam"
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}-apply.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}-apply.tsv"
    params:
        path=config["mainEnv"],
        do=r"{realign}",   # reAlign or noAlign
        extra="",  # optional
        java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["applyReAlign"]
    wrapper:
           config["wrapper"]+"/gatk3/indelrealigner"




rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads:config["threads"]["samtools"] 
    params: 
        path=config["mainEnv"],
        extra=""
    wrapper:
        config["wrapper"]+"samtools/index"


rule preBQSR: 
    input:
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.bam",
        bai="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.bam",
        ref=config["ref"]["genome"],
        known=[config["ref"]["1KGp1snp"],config["ref"]["dbsnp"],config["ref"]["1KGomni"],config["ref"]["hapmap"],config["ref"]["mills1KG"]],      
    output: 
        target_intervals="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR.intervals",
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_preBQSR.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_preBQSR.tsv"
    params:
        path=config["mainEnv"],
        extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ " ,  # optional
        java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["preBQSR"]
    wrapper:
           config["wrapper"]+"/gatk/baserecalibrator"

rule applyBQSR: 
    input:
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.bam",
        bai="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}.bam",
        ref=config["ref"]["genome"],
        target_intervals="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR.intervals"
    output: 
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR.bam",
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR.tsv"
    params:
        path=config["mainEnv"],
        extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ " ,  # optional
        java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["applyBQSR"]
    wrapper:
           config["wrapper"]+"/gatk/applyBQSR"
rule leftAlign:
    input:
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR.bam",
        bai="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR.bam.bai",
        ref=config["ref"]["genome"]
    output: 
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR_leftAlign.bam"
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR_leftAlign.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR_leftAlign.tsv"
    params:
        path=config["mainEnv"],
        extra=" --tmp-dir ../data/align/{aligner}/{case}/{sample}/ ",
        java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["leftAlign"]
    wrapper:
           config["wrapper"]+"/gatk/leftAlignIndels"
rule fixMate:
    input:
        ref=config["ref"]["genome"],
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR_leftAlign.bam",
        bai="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR_leftAlign.bam.bai"
    output: 
        bam="../data/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_BQSR_leftAlign_fixMate.bam"
    log:
        "../logs/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR_leftAlign_fixMate.logs"
    benchmark:
        "../benchmark/align/{aligner}/{case}/{sample}/{case}_{sample}_{qcpipe}_{aligner}_{rmdup}_{realign}_applyBQSR_leftAlign_fixMate.tsv"
    params:
        path=config["mainEnv"],
        extra=" --TMP_DIR ../data/align/{aligner}/{case}/{sample}/ ",
        java_opts=" -Xmx10G ",
    threads: config["threads"]["gatk"]["fixMate"]
    wrapper:
           config["wrapper"]+"/gatk/fixMateInformation"


