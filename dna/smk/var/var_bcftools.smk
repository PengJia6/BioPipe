# Samtools Calling
# https://samtools.github.io/bcftools/

rule Samtools_Mpileup:
    input:
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          vcf=path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.mpileup.g.vcf.gz"
    params:
          extra="",
          dp=5
    threads: config["threads"]["Samtools_Mpileup"]
    log:
       path_log + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.MergeVcf.tsv"
    run:
        shell(
            "{path_bcftools}bcftools mpileup -f {input.ref}  -g {params.dp} -O z -o {output.vcf} "
            "--threads {threads} {input.bam} 2>{log} 1>{log} ")

rule Samtools_Call:
    input:
         rules.Samtools_Mpileup.output.vcf
    output:
          vcf=path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.raw.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Samtools_Call"]
    log:
       path_log + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.Samtools_Call.logs"
    benchmark:
             path_bm + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.Samtools_Call.tsv"
    run:
        shell("{path_bcftools}bcftools call -mv -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

rule Samtools_Filter:
    input:
         rules.Samtools_Call.output.vcf
    output:
          vcf=path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Samtools_Filter"]
    log:
       path_log + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.Samtools_Filter.logs"
    benchmark:
             path_bm + "gremlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.Samtools_Filter.tsv"
    run:
        shell("{path_bcftools}bcftools view -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

# vcf=
