# Bcftools Calling
# https://samtools.github.io/bcftools/

rule Bcftools_Mpileup:
    input:
         unpack(getHQbamsample),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          vcf=path_data + "germlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.mpileup.g.vcf.gz"
    params:
          extra="",
          dp=5
    threads: config["threads"]["Bcftools_Mpileup"]
    log:
       path_log + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.MergeVcf.tsv"
    run:
        shell(
            "{path_bcftools}bcftools mpileup -f {input.ref}  -g {params.dp} -O z -o {output.vcf} "
            "--threads {threads} {input.bam} 2>{log} 1>{log} ")

rule Bcftools_Call:
    input:
         rules.Bcftools_Mpileup.output.vcf
    output:
          vcf=path_data + "germlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.raw.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Bcftools_Call"]
    log:
       path_log + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.Bcftools_Call.logs"
    benchmark:
             path_bm + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.Bcftools_Call.tsv"
    run:
        shell("{path_bcftools}bcftools call -mv -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

rule Bcftools_Filter:
    input:
         rules.Bcftools_Call.output.vcf
    output:
          vcf=path_data + "germlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["Bcftools_Filter"]
    log:
       path_log + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.Bcftools_Filter.logs"
    benchmark:
             path_bm + "gremlineVar/Bcftools/perSample/{bam_sample}/{bam_sample}.Bcftools.Bcftools_Filter.tsv"
    run:
        shell("{path_bcftools}bcftools view -i '%QUAL>=20' -Oz -o {output.vcf} {input} 2>{log} 1>{log} ")

