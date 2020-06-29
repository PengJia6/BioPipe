# def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
#     if regions:
#         params = "--intervals '{}' ".format(regions)
#         padding = config["processing"].get("region-padding")
#         if padding:
#             params += "--interval-padding {}".format(padding)
#         return params
#     return default
#
# def get_call_variants_params(wildcards, input):
#     return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
# #             config["params"]["gatk"]["HaplotypeCaller"])

rule HC_CallVar:
    input:
         unpack(getHQbamsample),
         # bai=unpack(),
         ref=path_genome,
         dict=path_dict,

    output:
          gvcf=path_data + "germlineVar/HC/perContig/{bam_sample}/{bam_sample}.{contig}.g.vcf.gz"
    log:
       path_log + "gremlineVar/HC/perContig/{bam_sample}/{bam_sample}.HC.{contig}.logs"
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
          path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.vcf.gz"
    params:
          extra="",
          java_options=""
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
    log:
       path_log + "gremlineVar/HC/jointCall/jointCall/joint.HC.MergeVcf.logs"
    benchmark:
             path_bm + "gremlineVar/HC/jointCall/jointCall/joint.HC.MergeVcf.tsv"
    run:
        # inputs = " ".join([("-V " + f) for f in input.gvcfs])
        inputs = " ".join(["INPUT={}".format(f) for f in input])
        shell("{path_picard}picard MergeVcfs {params.extra} "
              " {inputs} OUTPUT={output} 2>{log} 1>{log}")

#
# rule merge_variants:
#     input:
#         vcf=expand("../data/germlineVar/HC/jointCall/joint.{contig}.vcf.gz", contig=contigs)
#     output:
#         vcf="../data/germlineVar/HC/joint.HC.vcf.gz"
#     log:
#         "../logs/gremlineVar/HC/all.genotype.logs"
#     wrapper:
#         config["wrapper"]+"picard/mergevcfs"

"""
def get_vartype_arg(wildcards):
    return wildcards.vartype


rule VQSR:
    input:
        vcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.vcf.gz",
        ref=config["ref"]["genome"],
        # resources have to be given as named input files
        hapmap=config["ref"]["hapmap"],
        omni=config["ref"]["1KGomni"],
        g1k=config["ref"]["1KGp1snp"],
        dbsnp=config["ref"]["dbsnp"],
        mills=config["ref"]["mills1KG"],
        # use aux to e.g. download other necessary file
        aux=[config["ref"]["hapmap"]+".tbi",
		config["ref"]["1KGomni"]+".tbi",
        	config["ref"]["1KGp1snp"]+".tbi",
	        config["ref"]["dbsnp"]+".tbi" ,
		config["ref"]["mills1KG"]+".tbi",]
        
    output:
        vcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.{vartype}.recal.vcf",
        HQvcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.{vartype}.HQ.vcf",
        tranches="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.{vartype}.tranches",
	rscriptfile ="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.{vartype}.R"
    log:
        "../../../logs/vcf/gremline/HC/all.VQSR.{vartype}.logs"
    params:
        mode=get_vartype_arg,  # set mode, must be either SNP, INDEL or BOTH
        # resource parameter definition. Key must match named input files from above.
        resources={"hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
                   "omni":   {"known": False, "training": True, "truth": False, "prior": 12.0},
                   "g1k":   {"known": False, "training": True, "truth": False, "prior": 10.0},
                   "dbsnp":  {"known": True, "training": False, "truth": False, "prior": 2.0},
		   "mills":  {"known": False, "training": True, "truth": True, "prior": 14.0}
		},
        annotation=["QD", "FisherStrand","DP","FS","SOR","ReadPosRankSum","MQRankSum"],  # which fields to use with -an (see VariantRecalibrator docs)
        extra=" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 ",  # optional
        java_opts="", # optional
    wrapper:
        config["wrapper"]+"gatk/variantrecalibrator"

"""
