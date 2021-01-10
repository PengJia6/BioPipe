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
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.MergeVcf.tsv"
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
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC_JointCall_SelectSNV.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC_JointCall_SelectSNV.tsv"
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
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC_JointCall_SelectIndel.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
                 "name"] + ".HC_JointCall_SelectIndel.tsv"
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
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
           "name"] + ".HC_JointCall_FileterSNVHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
                 "name"] + ".HC_JointCall_FileterSNVHard.tsv"
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
         rules.HC_JointCall_SelectIndel.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.INDEL.passh.vcf.gz"
    threads: config["threads"]["HC_JointCall_FileterINDElHard"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
           "name"] + ".HC_JointCall_FileterINDElHard.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
                 "name"] + ".HC_JointCall_FileterINDElHard.tsv"
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
    input: rules.HC_MergeVCF.output
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
         rules.HC_SelectSNV.output
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
         rules.HC_SelectIndel.output
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

rule HC_Joint_MakeSitesOnlyVcf:
    input:
         rules.HC_MergeVCF_jointCall.output
    output:
          path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.sitesonly.vcf.gz"
    threads:
           config["threads"]["HC_Joint_MakeSitesOnlyVcf"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC_JointCall_sitesonly.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC_JointCall_sitesonly.tsv"
    run:
        shell("{path_gatk}gatk MakeSitesOnlyVcf -I {input}  -O {output} 2>{log} 1>{log}")

path_db_hapmap = config["ref"]["hapmap"]
path_db_1KGomni = config["ref"]["1KGomni"]
path_db_1KGp1snp = config["ref"]["1KGp1snp"]
path_db_dbsnp = config["ref"]["dbsnp"]
path_db_mills1KG = config["ref"]["mills1KG"]
path_db_axiomPoly = config["ref"]["axiomPoly"]
resources_index = [res + ".tbi" for res in
                   [path_db_dbsnp, path_db_mills1KG, path_db_1KGp1snp, path_db_hapmap, path_db_1KGomni,
                    path_db_axiomPoly]]


def get_vartype_paras(wildcards):
    if wildcards.vartype == "SNV":
        return "-mode SNP " \
               "-an QD -an MQRankSum  -an ReadPosRankSum -an FS -an SOR -an DP -an MQ " \
               "--max-gaussians 6"
    else:
        return "-mode INDEL " \
               "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP " \
               "--max-gaussians 4"


rule HC_Joint_VQSLOD:
    input:
         vcf=rules.HC_Joint_MakeSitesOnlyVcf.output,
         ref=path_genome,
         # resources have to be given as named input files
         hapmap=path_db_hapmap,
         omni=path_db_1KGomni,
         g1k=path_db_1KGp1snp,
         dbsnp=path_db_dbsnp,
         mills=path_db_mills1KG,
         axiomPoly=path_db_axiomPoly,
         resource_index=resources_index
    output:
          recal=path_data + "germlineVar/HC/jointCall/jointCall/" +
                config["project"]["name"] + ".HC.VQSLOD.{vartype}.recal",
          tranches=path_data + "germlineVar/HC/jointCall/jointCall/" +
                   config["project"]["name"] + ".HC.VQSLOD.{vartype}.tranches"
    threads:
           config["threads"]["HC_Joint_VQSLOD"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
           "name"] + ".HC_JointCall_Joint_VQSLOD_{vartype}.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
                 "name"] + ".HC_JointCall_Joint_VQSLOD_{vartype}.tsv"
    params:
          java_opt=" '-Xmx24g -Xms24g' ",
          vartype=get_vartype_paras
    run:
        if wildcards.vartype == "SNV":
            resources = \
                "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 " + input.hapmap + " " + \
                "--resource:omni,known=false,training=true,truth=false,prior=12.0 " + input.omni + " " + \
                "--resource:g1k,known=false,training=true,truth=false,prior=10.0 " + input.g1k + " " + \
                "--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 " + input.dbsnp + " "
        else:
            # if wildcards.vartype == "SNV":
            resources = \
                "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 " + input.dbsnp + " " + \
                "--resource:mills,known=false,training=true,truth=true,prior=12.0 " + input.mills + " " + \
                "--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 " + input.axiomPoly + " "
        shell("{path_gatk}gatk --java-options {params.java_opt}  VariantRecalibrator "
              "-V {input.vcf} "
              "--trust-all-polymorphic "
              "-R {input.ref} "
              "{params.vartype} "
              "-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 "
              "-tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 "
              "{resources} "
              "-O {output.recal} "
              "--tranches-file {output.tranches} "
              "2>{log} 1>{log}"
              )

rule HC_Joint_applyVQSR:
    input:
         vcf=rules.HC_MergeVCF_jointCall.output,
         recal=rules.HC_Joint_VQSLOD.output.recal,
         tranches=rules.HC_Joint_VQSLOD.output.tranches,
    output:
          vcf=path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"][
              "name"] + ".HC.{vartype}.VQSR.vcf.gz"

    threads:
           config["threads"]["HC_Joint_applyVQSR"]
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
           "name"] + ".HC_JointCall_Joint_applyVQSR_{vartype}.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/" + config["project"][
                 "name"] + ".HC_JointCall_Joint_applyVQSR_{vartype}.tsv"
    params:
          java_opt=" '-Xmx5g -Xms5g' ",
    run:
        vartype = "SNP" if wildcards.vartype == "SNV" else "INDEL"
        shell("{path_gatk}gatk --java-options {params.java_opt}  ApplyVQSR "
              "-V {input.vcf}  --recal-file {input.recal} --tranches-file {input.tranches} "
              "--truth-sensitivity-filter-level 99.7 "
              # "--create-output-variant-index fasle "
              "-mode {vartype} "
              "-O {output.vcf} "
              "2>{log} 1>{log}")

rule HC_VQSLOD:
    input:
         vcf=rules.HC_MergeVCF.output,
         ref=path_genome,
         # resources have to be given as named input files
         hapmap=path_db_hapmap,
         omni=path_db_1KGomni,
         g1k=path_db_1KGp1snp,
         dbsnp=path_db_dbsnp,
         mills=path_db_mills1KG,
         axiomPoly=path_db_axiomPoly,
         resource_index=resources_index
    output:
          recal=path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.VQSLOD.{vartype}.recal",
          tranches=path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.VQSLOD.{vartype}.tranches"

    threads:
           config["threads"]["HC_Joint_VQSLOD"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_JointCall_Joint_VQSLOD_{vartype}.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jperSample/{bam_sample}/{bam_sample}.HC_JointCall_Joint_VQSLOD_{vartype}.tsv"
    params:
          java_opt=" '-Xmx10g -Xms10g' ",
          vartype=get_vartype_paras
    run:
        if wildcards.vartype == "SNV":
            resources = \
                "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 " + input.hapmap + " " + \
                "--resource:omni,known=false,training=true,truth=false,prior=12.0 " + input.omni + " " + \
                "--resource:g1k,known=false,training=true,truth=false,prior=10.0 " + input.g1k + " " + \
                "--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 " + input.dbsnp + " "
        else:
            # if wildcards.vartype == "SNV":
            resources = \
                "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 " + input.dbsnp + " " + \
                "--resource:mills,known=false,training=true,truth=true,prior=12.0 " + input.mills + " " + \
                "--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 " + input.axiomPoly + " "
        shell("{path_gatk}gatk --java-options {params.java_opt}  VariantRecalibrator "
              "-V {input.vcf} "
              "--trust-all-polymorphic "
              "-R {input.ref} "
              "{params.vartype} "
              "-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 "
              "-tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 "
              "{resources} "
              "-O {output.recal} "
              "--tranches-file {output.tranches} "
              "2>{log} 1>{log}"
              )

rule HC_applyVQSR:
    input:
         vcf=rules.HC_MergeVCF.output,
         recal=rules.HC_VQSLOD.output.recal,
         tranches=rules.HC_VQSLOD.output.tranches,
    output:
          vcf=path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.{vartype}.VQSR.vcf.gz"

    threads:
           config["threads"]["HC_applyVQSR"]
    log:
       path_log + "gremlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC_JointCall_Joint_applyVQSR_{vartype}.logs"
    benchmark:
             path_bm + "gremlineVar/HC/j/perSample/{bam_sample}/{bam_sample}.HC_JointCall_Joint_applyVQSR_{vartype}.tsv"
    params:
          java_opt=" '-Xmx5g -Xms5g' ",
    run:
        vartype = "SNP" if wildcards.vartype == "SNV" else "INDEL"
        shell("{path_gatk}gatk --java-options {params.java_opt}  ApplyVQSR "
              "-V {input.vcf}  --recal-file {input.recal} --tranches-file {input.tranches} "
              "--truth-sensitivity-filter-level 99.7 "
              # "--create-output-variant-index fasle "
              "-mode {vartype} "
              "-O {output.vcf} "
              "2>{log} 1>{log}")
