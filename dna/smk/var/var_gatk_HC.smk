rule HC_CallVar:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome+".fai",
         sindex2=path_genome+".fai_chrom",
         dict=path_dict,


    output:
          gvcf=path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{contig}.g.vcf.gz"
    log:
       path_log + "gremlineVar/HC/perContig/{bam_sample}/{bam_sample}.HC.{contig}.logs"
    threads: config["threads"]["HC_CallVar"]
    benchmark:
             path_bm + "gremlineVar/HC/perContig/{bam_sample}/{bam_sample}.HC.{contig}.tsv"
    params:
          extra="",
          java_options="",
          regions="",
          dbsnp=[],
    run:
        thispara = "" if len(params.dbsnp) == 0 else " --dbsnp {params.dbsnp[0]}"
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{path_gatk}gatk {java_opt} HaplotypeCaller {params.extra} {thispara} "
              " -R {input.ref} -ERC GVCF -L {wildcards.contig} -I {input.bam} -O {output.gvcf}"
              " 2>{log} 1>{log}")

rule HC_GT:
    input:
         ref=path_genome,
         gvcf=path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{contig}.g.vcf.gz"
    output:
          vcf=path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{contig}.vcf.gz"
    params:
          extra="",
          java_options=""
          # extra=config["params"]["gatk"]["GenotypeGVCFs"]
    threads: config["threads"]["HC_GT"]

    log:
       path_log + "gremlineVar/HC/perContig/{bam_sample}/{bam_sample}.HC.genotype_variants.{contig}.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perContig/{bam_sample}/{bam_sample}.HC.genotype_variants.{contig}.tsv"
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "

        shell("{path_gatk}gatk {java_opt}  GenotypeGVCFs {params.extra} "
              " -R {input.ref} -L {wildcards.contig} -V {input.gvcf} -O {output.vcf} 2>{log} 1>{log}")
rule HC_MergeVCF:
    input:
         expand(path_data + "germlineVar/HC/perContig/{{bam_sample}}/{{bam_sample}}.{contig}.vcf.gz", contig=contigs)
    output:
          path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.raw.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["HC_MergeVCF"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.MergeVcf.tsv"
    run:
        inputs = " ".join(["INPUT={}".format(f) for f in input])
        shell("{path_picard}picard MergeVcfs {params.extra} "
              " {inputs} OUTPUT={output} 2>{log} 1>{log}")

rule HC_CombineCalls:
    input:
         ref=path_genome,
         gvcfs=expand(path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{{contig}}.g.vcf.gz",
                      bam_sample=bam_sample_list)
    output:
          path_data + "germlineVar/HC/jointCall/perContig/{contig}.g.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["HC_CombineCalls"]
    log:
       path_log + "gremlineVar/HC/jointCall/perContig/{contig}.CombineCalls.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/perContig/{contig}.CombineCalls.tsv"
    run:
        # inputs = " ".join(expand(["-V {vcf}"],vcf={input.gvcfs}))
        inputs = " ".join([("-V " + f) for f in input.gvcfs])

        # inputs=" ".join([ str("-V "+str(i)) for i in {input.gvcfs}])
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "

        shell("{path_gatk}gatk {java_opt} CombineGVCFs {params.extra} "
              " {inputs} -R {input.ref} -O {output} 2>{log} 1>{log}")

rule HC_GT_jointCall:
    input:
         ref=path_genome,
         gvcf=path_data + "germlineVar/HC/jointCall/perContig/{contig}.g.vcf.gz"
    output:
          vcf=path_data + "germlineVar/HC/jointCall/perContig/{contig}.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["HC_GT_jointCall"]
    log:
       path_log + "gremlineVar/HC/jointCall/perContig/{contig}.GT.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/perContig/{contig}.GT.tsv"
    run:
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{path_gatk}gatk {java_opt}  GenotypeGVCFs {params.extra} "
              " -R {input.ref} -V {input.gvcf} -O {output.vcf} 2>{log} 1>{log}")

rule HC_MergeVCF_jointCall:
    input:
         expand(path_data + "germlineVar/HC/jointCall/perContig/{contig}.vcf.gz", contig=contigs)
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.vcf.gz"
    params:
          extra="",
          java_options=""
    threads: config["threads"]["HC_MergeVCF_jointCall"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC.MergeVcf.tsv"
    run:
        # inputs = " ".join([("-V " + f) for f in input.gvcfs])
        inputs = " ".join(["INPUT={}".format(f) for f in input])
        shell("{path_picard}picard MergeVcfs {params.extra} "
              " {inputs} OUTPUT={output} 2>{log} 1>{log}")

