localrules: all
configfile: "conf/config.yaml"

import pandas as pd 
contigs = pd.read_table(config["ref"]["genome"] + ".fai",
		     	header=None,usecols=[0] , squeeze = True, dtype = str)
contigs = contigs.loc[contigs.index[:3]]

samples = pd.read_csv(config["bam"], index_col=0, dtype=str)

def get_one_sample_bam(wildcards):
    sample=wildcards.sample
    bam=samples.loc[sample,"bampath"]
    bai=bam+".bai"
    return {"bam":bam, "bai":bai}

wildcard_constraints:
    vartype="SNP|INDEL",
    contig="|".join(contigs),
    sample="|".join(samples.index)



rule all: 
    input: 
       [ "../data/HQBam/{sample}.bam".format(sample = sample) for sample in samples.index ],
       [ "../data/HQBam/{sample}.bam.bai".format(sample = sample) for sample in samples.index ],
       "../data/germlineVar/HC/joint.HC.vcf.gz"

rule loadbam:
    input:
       unpack(get_one_sample_bam)
    
    output:
       bam="../data/HQBam/{sample}.bam",
       bai="../data/HQBam/{sample}.bam.bai"
    shell: 
       """
       ln -sr {input.bam} {output.bam} 
       ln -sr {input.bai} {output.bai}
       """


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
        bam="../data/HQBam/{sample}.bam",
        bai="../data/HQBam/{sample}.bam.bai",
        ref=config["ref"]["genome"],
        regions=[]
    output:
        gvcf="../data/germlineVar/HC/jointCall/{sample}/{sample}.{contig}.g.vcf.gz"
    log:
        "../logs/gremlineVar/HC/jointCall/{sample}/{sample}.HC.{contig}.logs"
    params:
        extra=get_call_variants_params
    wrapper: config["wrapper"]+"gatk/haplotypecaller"

rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("../data/germlineVar/HC/jointCall/{sample}/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="../data/germlineVar/HC/jointCall/joint.{contig}.g.vcf.gz"
    log:
        
        "../logs/gremlineVar/HC/joint.HC.{contig}.logs"
    wrapper:
        config["wrapper"]+"gatk/combinegvcfs"

rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="../data/germlineVar/HC/jointCall/joint.{contig}.g.vcf.gz"
    output:
        vcf="../data/germlineVar/HC/jointCall/joint.{contig}.vcf.gz"
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "../logs/gremlineVar/HC/joint.genotype.{contig}.logs"
    wrapper:
        config["wrapper"]+"gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcf=expand("../data/germlineVar/HC/jointCall/joint.{contig}.vcf.gz", contig=contigs)
    output:
        vcf="../data/germlineVar/HC/joint.HC.vcf.gz"
    log:
        "../logs/gremlineVar/HC/all.genotype.logs"
    wrapper:
        config["wrapper"]+"picard/mergevcfs"

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
