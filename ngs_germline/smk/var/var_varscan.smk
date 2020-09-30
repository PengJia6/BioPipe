# Samtools Calling
# https://samtools.github.io/bcftools/

rule Samtools_Mpileup:
    input:
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.Samtools.mpileup"
    params:
          extra="",
          dp=5
    threads: config["threads"]["Samtools_Mpileup"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Samtools.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Samtools.MergeVcf.tsv"
    run:
        shell(
            "{path_samtools}samtools mpileup -B -f {input.ref} -o {output} {input.bam} "
            "2>{log} 1>{log} ")

rule Varscan_Call:
    input:
         rules.Samtools_Mpileup.output
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.raw.vcf.gz"
    params:
          extra=" --p-value 0.05 ",
    threads: config["threads"]["Varscan_Call"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call.tsv"
    run:
        shell("{path_varscan}varscan mpileup2cns {input} {params.extra}  --output-vcf 1 | "
              "{path_bcftools}bcftools view -Oz -o {output} 1>{log} 2>{log} ")

rule Varscan_Filter:
    input:
         rules.Varscan_Call.output.vcf
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Varscan_Filter"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter.tsv"
    run:
        shell("{path_bcftools}bcftools view {params.extra} -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

rule Varscan_Call_SNV:
    input:
         rules.Samtools_Mpileup.output
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.SNV.raw.vcf.gz"
    params:
          extra=" --p-value 0.05 ",
    threads: config["threads"]["Varscan_Call_SNV"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call_SNV.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call_SNV.tsv"
    run:
        shell("{path_varscan}varscan mpileup2snp {input} {params.extra} --output-vcf 1 | "
              "{path_bcftools}bcftools view -Oz -o {output} 1>{log} 2>{log} ")

rule Varscan_Filter_SNV:
    input:
         rules.Varscan_Call_SNV.output.vcf
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.SNV.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Varscan_Filter_SNV"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter_SNV.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter_SNV.tsv"
    run:
        shell("{path_bcftools}bcftools view {params.extra} -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

rule Varscan_Call_INDEL:
    input:
         rules.Samtools_Mpileup.output
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.INDEL.raw.vcf.gz"
    params:
          extra=" --p-value 0.05 ",
    threads: config["threads"]["Varscan_Call_INDEL"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call_INDEL.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Call_INDEL.tsv"
    run:
        shell("{path_varscan}varscan mpileup2indel {input}  {params.extra} --output-vcf 1 | "
              "{path_bcftools}bcftools view -Oz -o {output} 1>{log} 2>{log} ")

rule Varscan_Filter_INDEL:
    input:
         rules.Varscan_Call_INDEL.output.vcf
    output:
          vcf=path_data + "germlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.INDEL.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Varscan_Filter_INDEL"]
    log:
       path_log + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter_INDEL.logs"
    benchmark:
             path_bm + "gremlineVar/varscan/perSample/{bam_sample}/{bam_sample}.varscan.Varscan_Filter_INDEL.tsv"
    run:
        shell("{path_bcftools}bcftools view {params.extra} -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")
