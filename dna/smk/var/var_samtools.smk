# Samtools Calling
# https://samtools.github.io/bcftools/
# http://arxiv.org/abs/1207.3907

rule Samtools_Mpileup:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          mpileup=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf",
          vcfgz=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf.gz"
    params:
          extra="",
          chunks=100000000
    threads: config["threads"]["FB_CallVar"]
    log:
       path_log + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.MergeVcf.tsv"
    run:
        shell("{path_freebayes}freebayes-parallel <({path_freebayes}fasta_generate_regions.py "
              "{input.ref} {params.chunks}) {threads} "
              " -f {input.ref} {input.bam} > {output.vcf} 2>{log}")
        shell("{path_bgzip}bgzip -c {output.vcf} > {output.vcfgz}")
#
rule FB_Filter:
    input:
         vcf=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf"
    output:
          vcf=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf",
          vcfgz=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["FB_Filter"]
    log:
       path_log + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.FB_Filter.logs"

    benchmark:
             path_bm + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.FB_Filter.tsv"
    run:
        shell("{path_vcflib}vcffilter -f 'QUAL > 20' {input.vcf}>{output.vcf}")
        shell("{path_bgzip}bgzip -c {output.vcf} > {output.vcfgz}")
