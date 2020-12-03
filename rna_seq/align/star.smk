# ======================================================================================================================
# Project: BioPipe
# Script : star.smk TODO check 
# Author : Peng Jia
# Date   : 2020.11.29
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
rule STAR:
    input:
         R1=path_data + "cleandata/{qc_pipe}/{sample}/{sample}_{unit}_{qc_pipe}_1.fq.gz",
         R2=path_data + "cleandata/{qc_pipe}/{sample}/{sample}_{unit}_{qc_pipe}_2.fq.gz",
         gtf=path_gtf
    output:
          tag=path_data + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR",
          bam=path_data + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR_Aligned.sortedByCoord.out.bam",
          junction=path_data + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR_Chimeric.out.junction",

    log:
       path_log + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR.log"
    benchmark:
             path_log + "align/STAR/{sample}/{sample}_{unit}_{qc_pipe}_STAR.tsv"
    threads: config["threads"]["STAR"]
    params:
          extra="",
          read_len=get_read_len
    run:
        param_read_len = int(params.read_len) - 1
        shell("{path_star} --runThreadN {threads} "
              "--genomeDir {path_STAR_index}_{param_read_len} --readFilesIn {input.R1} {input.R2} "
              "--readFilesCommand zcat --sjdbGTFfile {input.gtf} --sjdbOverhang {param_read_len} "
              "--outFileNamePrefix {output.tag}_ --outSAMtype BAM SortedByCoordinate "
              "--outSAMstrandField intronMotif "
              "--outSAMunmapped Within "
              "--chimSegmentMin 12 "  # ** essential to invoke chimeric read detection & reporting **
              "--chimJunctionOverhangMin 8 "
              "--chimOutJunctionFormat 1 "  # **essential** includes required metadata in Chimeric.junction.out file.
              "--alignSJDBoverhangMin 10 "
              "--alignMatesGapMax 100000 "  # avoid readthru fusions within 100k
              "--alignIntronMax 100000 "
              "--alignSJstitchMismatchNmax 5 -1 5 5 "  # settings improved certain chimera detections
              "--outSAMattrRGline ID:GRPundef "
              "--chimMultimapScoreRange 3 "
              "--chimScoreJunctionNonGTAG -4 "
              "--chimMultimapNmax 20 "
              "--chimNonchimScoreDropMin 10 "
              "--peOverlapNbasesMin 12 "
              "--peOverlapMMp 0.1 "
              "--alignInsertionFlush Right "
              "--alignSplicedMateMapLminOverLmate 0 "
              "--alignSplicedMateMapLmin 30 "
              "--outBAMsortingThreadN {threads} --quantMode TranscriptomeSAM GeneCounts 2>{log} 1>{log}")
        shell("touch {output.tag}")

rule MergeBam:
    input:
         get_sample_bams
         #        expand("{path_prefix}align/{{aligner}}/{{case}}/{{sample}}/{{case}}_{{sample}}_{unit}_{{qcpipe}}_{{aligner}}_sorted.bam", unit=caseinfo.loc[({case},{sample}),"unit"],path_prefix=[path_data])
    output:
          bam=path_data + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}.bam",
    log:
       path_log + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}.merge.logs"
    benchmark:
             path_log + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}.merge.tsv"
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

rule MergeJunction:
    input:
         get_sample_jucntion
         #expand("{path_prefix}align/{{aligner}}/{{case}}/{{sample}}/{{case}}_{{sample}}_{unit}_{{qcpipe}}_{{aligner}}_sorted.bam", unit=caseinfo.loc[({case},{sample}),"unit"],path_prefix=[path_data])
    output:
          junction=path_data + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}_Chimeric.out.junction"
    log:
       path_log + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}.merge_junction.logs"
    benchmark:
             path_log + "align/{aligner}/{sample}/{sample}_{qcpipe}_{aligner}.merge_junction.tsv"
    threads: config["threads"]["MergeJunction"]
    params:
          # path=path_samtools,
          extra=" ",
    run:
        if len(input) < 2:
            shell("ln -sr {input} {output.junction} ")
            shell("sleep 1")
            shell("touch -h {output.junction} 1>>{log} 2>>{log}")
            shell("echo only one bam, make soft link 1>>{log} 2>>{log}")
        else:
            file = open(output.junction, "w")
            num = 0
            for item in output:
                num += 1
                if num == 1:
                    for line in open(item):
                        file.write(line)
                else:
                    line_num = 0
                    for line in open(item):
                        line_num += 1
                        if line_num == 1: continue
                        file.write(line)
            file.close()

            # shell("{path_samtools}samtools merge -@ {threads} {params.extra} "
            #       "{output.junction} {input} 2>{log} 1>{log}")
