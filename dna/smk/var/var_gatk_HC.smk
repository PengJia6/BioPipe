# Gatk HaplotypeCaller
# https://gatk.broadinstitute.org/
rule HC_CallVar:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         sindex=path_genome + ".fai",
         sindex2=path_genome + ".fai_chrom",
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
         # gvcf=path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{contig}.g.vcf.gz"
         gvcf=rules.HC_CallVar.output.gvcf
    output:
          # vcf=rules.HC_CallVar.output.gvcf
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
         # gvcf=path_data + "germlineVar/HC/jointCall/perContig/{contig}.g.vcf.gz"
         gvcf=rules.HC_CombineCalls.output
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
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.raw.vcf.gz"
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

## Filter variant
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

rule HC_JointCall_SelectSNV:
    input:
         rules.HC_MergeVCF_jointCall.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.SNV.vcf.gz"
    threads:
           config["threads"]["HC_JointCall_SelectSNV"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_SelectSNV.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_SelectSNV.tsv"
    run:
        shell("{path_gatk}gatk SelectVariants -V {input} -select-type SNP -O {output}"
              " 2>{log} 1>{log}")
rule HC_JointCall_SelectIndel:
    input:
         rules.HC_MergeVCF_jointCall.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.INDEL.vcf.gz"
    threads:
           config["threads"]["HC_JointCall_SelectIndel"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_SelectIndel.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_SelectIndel.tsv"
    run:
        shell("{path_gatk}gatk SelectVariants -V {input} -select-type INDEL -O {output}"
              " 2>{log} 1>{log}")

rule HC_JointCall_FileterSNVHard:
    input:
         rules.HC_JointCall_SelectSNV.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.SNV.passh.vcf.gz"
    threads: config["threads"]["HC_JointCall_FileterSNVHard"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_FileterSNVHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_FileterSNVHard.tsv"
    run:
        shell('{path_gatk}gatk VariantFiltration '
              ' -V {input} '
              ' -filter "QD < 2.0"  --filter-name "QD2" '
              ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
              ' -filter "SOR > 3.0" --filter-name "SOR3" '
              ' -filter "FS > 60.0" --filter-name "FS60" '
              ' -filter "MQ < 40.0" --filter-name "MQ40" '
              ' -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
              ' -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"  '
              ' -O {output} '
              ' 2>{log} 1>{log}')

rule HC_JointCall_FileterINDElHard:
    input:
         rules.HC_JointCall_SelectSNV.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.INDEL.passh.vcf.gz"
    threads: config["threads"]["HC_JointCall_FileterINDElHard"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_FileterINDElHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC_JointCall_FileterINDElHard.tsv"
    run:
        shell('{path_gatk}gatk VariantFiltration '
              ' -V {input} '
              ' -filter "QD < 2.0"  --filter-name "QD2" '
              ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
              ' -filter "FS > 200.0" --filter-name "FS200" '
              ' -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"  '
              ' -O {output} '
              ' 2>{log} 1>{log}')

########################################################################################################################
####### Hard Filtering for single sample
rule HC_SelectSNV:
    input: rules.HC_MergeVCF.output
    output:
          path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.SNV.vcf.gz"
    threads: config["threads"]["HC_SelectSNV"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_SelectSNV.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_SelectSNV.tsv"
    run:
        shell("{path_gatk}gatk SelectVariants -V {input} -select-type SNP -O {output}"
              " 2>{log} 1>{log}")
rule HC_SelectIndel:
    input: rules.HC_MergeVCF_jointCall.output
    output:path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.INDEL.vcf.gz"

    threads: config["threads"]["HC_SelectIndel"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_SelectIndel.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_SelectIndel.tsv"
    run:
        shell("{path_gatk}gatk SelectVariants -V {input} -select-type INDEL -O {output}"
              " 2>{log} 1>{log}")

rule HC_FileterSNVHard:
    input:
         rules.HC_JointCall_SelectSNV.output
    output:
          path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.SNV.passh.vcf.gz"
    threads: config["threads"]["HC_FileterSNVHard"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterSNVHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterSNVHard.tsv"
    run:
        shell('{path_gatk}gatk VariantFiltration '
              ' -V {input} '
              ' -filter "QD < 2.0"  --filter-name "QD2" '
              ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
              ' -filter "SOR > 3.0" --filter-name "SOR3" '
              ' -filter "FS > 60.0" --filter-name "FS60" '
              ' -filter "MQ < 40.0" --filter-name "MQ40" '
              ' -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
              ' -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"  '
              ' -O {output} '
              ' 2>{log} 1>{log}')

rule HC_FileterINDElHard:
    input:
         rules.HC_JointCall_SelectSNV.output
    output:
          path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.INDEL.passh.vcf.gz"

    threads: config["threads"]["HC_FileterINDElHard"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterINDElHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterINDElHard.tsv"
    run:
        shell('{path_gatk}gatk VariantFiltration '
              ' -V {input} '
              ' -filter "QD < 2.0"  --filter-name "QD2" '
              ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
              ' -filter "FS > 200.0" --filter-name "FS200" '
              ' -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"  '
              ' -O {output} '
              ' 2>{log} 1>{log}')
# rule HC_FileterINDElHard:
#     input:
#          rules.HC_JointCall_SelectSNV.output
#     output:
#           path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.INDEL.passh.vcf.gz"
#
#     threads: config["threads"]["HC_FileterINDElHard"]
#     log:
#        path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterINDElHard.logs"
#     benchmark:
#              path_bm + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_FileterINDElHard.tsv"
#     run:
#         shell('{path_gatk}gatk VariantFiltration '
#               ' -V {input} '
#               ' -filter "QD < 2.0"  --filter-name "QD2" '
#               ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
#               ' -filter "FS > 200.0" --filter-name "FS200" '
#               ' -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"  '
#               ' -O {output} '
#               ' 2>{log} 1>{log}')

