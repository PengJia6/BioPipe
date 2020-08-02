
localrules: all
pindel_types = ["D", "BP", "INV", "TD", "LI", "SI", "RP"]
configfile: "case.yaml"
samples= ["CC","CN","CP","CPC","LM","PBMC","RC","RN"]
import pandas as pd 
contigs= pd.read_table(config["ref"]["genome"]+".fai2",
		     	header=None,usecols=[0],squeeze=True,dtype=str)


rule all: 
    input: 
#       ["../../../data/bam/bysample/{sample}/{sample}.rmdup.indelrealign.BQSR.bam.bai".format(sample=sample) for sample in config["all"]  ]+
 #      ["../../../data/vcf/germline/{sample}/{sample}.rmdup.indelrealign.BQSR.{contig}.g.vcf.gz".format(sample=sample,contig=contig) for sample in config["all"] for contig in contigs], 
        "../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.vcf.gz",
        "../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.SNP.HQ.vcf",
        "../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.VQSR.INDEL.HQ.vcf"





wildcard_constraints:
    vartype="SNP|INDEL",
    contig="|".join(contigs)

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])

rule call_variants:
    input:
        bam="../../../data/bam/bysample/{sample}/{sample}.rmdup.indelrealign.BQSR.bam",
        ref=config["ref"]["genome"],
        known=config["ref"]["dbsnp"],
        regions=[]
    output:
        gvcf=protected("../../../data/vcf/germline/HC/{sample}/{sample}.rmdup.indelrealign.BQSR.{contig}.g.vcf.gz")
    log:
        "../../../logs/vcf/gremline/HC/{sample}/{sample}.HC.{contig}.logs"
    params:
        extra=get_call_variants_params
    wrapper: config["wrapper"]+"gatk/haplotypecaller"

rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("../../../data/vcf/germline/HC/{sample}/{sample}.rmdup.indelrealign.BQSR.{{contig}}.g.vcf.gz", sample=samples)
    output:
        gvcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.{contig}.g.vcf.gz"
    log:
        
        "../../../logs/vcf/gremline/HC/all.HC.{contig}.logs"
    wrapper:
        config["wrapper"]+"gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.{contig}.g.vcf.gz"
    output:
        vcf=temp("../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "../../../logs/vcf/gremline/HC/all.genotype.{contig}.logs"
    wrapper:
        config["wrapper"]+"gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcf=expand("../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.{contig}.vcf.gz", contig=contigs)
    output:
        vcf="../../../data/vcf/germline/HC/all.rmdup.indelrealign.BQSR.HC.vcf.gz"
    log:
        "../../../logs/vcf/gremline/HC/all.genotype.logs"
    wrapper:
        config["wrapper"]+"picard/mergevcfs"


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


