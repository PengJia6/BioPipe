# Freebayes Calling
# https://github.com/ekg/freebayes
# http://arxiv.org/abs/1207.3907
rule FB_CallVar:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai"
    output:
          vcf=temp(path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf"),
          vcfgz=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf.gz"
    params:
          extra="",
          chunks=100000000
    threads: config["threads"]["FB_CallVar"]
    log:
       path_log + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.MergeVcf.tsv"
    shell:
         ## [Warning] freebayes-parallel need to PATH env
         """
         export PATH={path_freebayes}:$PATH
         {path_freebayes}freebayes-parallel <({path_freebayes}fasta_generate_regions.py {input.ref} {params.chunks}) {threads} -f {input.ref} {input.bam} > {output.vcf} 2>{log} 
         {path_bgzip}bgzip -c {output.vcf} > {output.vcfgz} 
         """
#         shell("{path_freebayes}freebayes-parallel <({path_freebayes}fasta_generate_regions.py "
#       "{input.ref} {params.chunks}) {threads} "
#       " -f {input.ref} {input.bam} > {output.vcf} 2>{log}")
# shell("{path_bgzip}bgzip -c {output.vcf} > {output.vcfgz}")
#
rule FB_Filter:
    input:
         rules.FB_CallVar.output.vcfgz
         # vcf=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf"
    output:
          # vcf=temp(path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf"),
          vcfgz=path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf.gz"
    params:
          extra="",
    threads: config["threads"]["FB_Filter"]
    log:
       path_log + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.FB_Filter.logs"

    benchmark:
             path_bm + "gremlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.FB_Filter.tsv"
    run:
        shell("{path_bcftools}bcftools view -i '%QUAL>=20' -Oz -o {output.vcfgz} {input} 2>{log} 1>{log} ")

        # shell("{path_vcflib}vcffilter -f 'QUAL > 20' {input.vcf}>{output.vcf}")
        # shell("{path_bgzip}bgzip -c {output.vcf} > {output.vcfgz}")
