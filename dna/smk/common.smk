localrules: LoadGenome
from snakemake.io import expand
import pandas as pd
import os
from snakemake.shell import shell

# report: "../reports/workflow.rst"
configfile: "conf/config.yaml"
path_data = config["path_data"]
path_data = path_data if path_data[-1] == "/" else path_data + "/"
path_log = config["path_log"]
path_log = path_log if path_log[-1] == "/" else path_log + "/"
path_bm = config["path_bm"]
path_bm = path_bm if path_bm[-1] == "/" else path_bm + "/"

# validate(config, schema="../schemas/config.schema.yaml")
caseinfo = pd.read_csv(config["caseinfo"]).set_index(["case", "sample", "unit"], drop=False)
# sampleinfo= pd.read_csv(config["caseinfo"]).set_index(["case","sample"], drop=False)
# validate(samples, schema="../schemas/samples.schema.yaml")

# units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
caseinfo.index = caseinfo.index.set_levels([i.astype(str) for i in caseinfo.index.levels])  # enforce str in index

# validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai_chrom",
                        header=None, usecols=[0], squeeze=True, dtype=str)
#
# casedict = {}
# for index, info in caseinfo.iterrows():
#     case = index[0]
#     sample = index[1]
#     if case not in casedict:
#         casedict[case] = {}
#     casedict[case][sample] = ""
# num = 0
# sampleinfo = pd.DataFrame()
# for case in casedict:
#     for sample in casedict[case]:
#         num += 1
#         sampleinfo.loc[num, "case"] = case
#         sampleinfo.loc[num, "sample"] = sample
#         sampleinfo.loc[num, "id"] = case + sample

# sampleinfo = sampleinfo.set_index(["case", "sample"])
# sampleinfo.index = sampleinfo.index.set_levels([i.astype(str) for i in sampleinfo.index.levels])  # enforce str in index
# for case in casedict:
#     for sample in casedict[case]:
#         num += 1
#         sampleinfo.loc[(case, sample), "case"] = case
#         sampleinfo.loc[(case, sample), "sample"] = sample

bam_suffix = "_".join([config["pipe"]["trim"], config["pipe"]["aligner"], config["pipe"]["markdup"],
                       config["pipe"]["realign"], config["pipe"]["bqsr"], config["pipe"]["leftAlign"],
                       config["pipe"]["fixMate"], ])
path_HQbamList = expand([path_data + "HQbam/{u.case}_{u.sample}.bam.bai"], u=caseinfo.itertuples())
bam_sample_list = []
sampleinfo = pd.DataFrame()
if config["loadBamFormat"] == "csv":
    sampleinfo = pd.read_csv(config["bamsample"], index_col=0)
    bam_sample_list = list(sampleinfo.index)
else:
    bam_sample_list = list(set(expand(["{u.case}_{u.sample}"], u=caseinfo.itertuples())))
    for one in bam_sample_list:
        sampleinfo.loc[one, "path"] = path_data + "HQbam/{bamsample}.bam.bai".format(bamsample=one)


def getHQbamsample(wildcards):
    if config["loadBamFormat"] == "csv":
        sampleinfo = pd.read_csv(config["bamsample"], index_col=0)
        bam = sampleinfo.loc[wildcards.bam_sample, "path"]
        return {"bam": bam, "bai": bam + ".bai"}

    else:
        bam = path_data + "HQbam/{bam_sample}.bam".format(bam_sample=wildcards.bam_sample)

        return {"bam": bam, "bai": bam + ".bai"}


path_raw_vcf = []
if "HC" in config["pipe"]["snpindel"]:
    # joint calling
    path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.vcf.gz")
    # single sample calling
    for i in bam_sample_list:
        path_raw_vcf.append(
            path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.vcf.gz".format(bam_sample=i))
##### Wildcard constraints #####
wildcard_constraints:
                    vartype="snvs|indels",
                    qcpipe="fastp|trim|passqc",
                    aligner="bwa|bowtie",
                    markdup="bioMarkDup|picardMarkDup|noneMarkDup",
                    realign="reAlign|noneReAlign",
                    bqsr="noneBQSR|BQSR",
                    leftAlign="noneLeftAlign|LeftAlign",
                    fixMate="noneFixMate|FixMate",
                    R="1|2",
                    case="|".join(caseinfo["case"]),
                    sample="|".join(caseinfo["sample"]),
                    unit="|".join(caseinfo["unit"]),
                    LB="|".join(caseinfo["LB"]),
                    contig="|".join(contigs),
                    bam_sample="|".join(bam_sample_list)

                    # fuction for input and params


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"R1": fastqs.fq1, "R2": fastqs.fq2}
    return {"R1": fastqs.fq1, "R2": fastqs.fq2}


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{case}_{sample}\tSM:{case}_{sample}\tPL:{platform}\tLB:{LB}'".format(
        sample=wildcards.sample,
        case=wildcards.case,
        platform=caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), "PL"][0],
        LB=caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), "LB"][0])


# for bam_merge input
def get_sample_bam(wildcards):
    units = list(set(caseinfo.loc[(wildcards.case, wildcards.sample), "unit"]))
    aligner = wildcards.aligner
    case = wildcards.case
    sample = wildcards.sample
    qcpipe = wildcards.qcpipe
    return [path_data + "align/" + case + "/" + sample + "/" + case + "_" + sample + "_" + str(
        i) + "_" + qcpipe + "_" + aligner + "_sorted.bam" for i in units]


# return expand(["../data/align/"+wildcards.aligner+"/"+wildcards.case+"/"+wildcards.sample+"/"+wildcards.case+"_"+wildcards.sample+"_{unitD}_"+wildcards.qcpipe+"_"+wildcards.aligner+"_sorted.bam" ],unitD=["1","2"])
# return expand(["../data/align/"+wildcards.aligner+"/"+wildcards.case+"/"+wildcards.sample+"/"+wildcards.case+"_"+wildcards.sample+"_{unitD}_"+wildcards.qcpipe+"_"+wildcards.aligner+"_sorted.bam" ],unitD=caseinfo.loc[(wildcards.case,wildcards.sample)].unit)


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams2222(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)


#
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
#             config["params"]["gatk"]["HaplotypeCaller"])
#
#
# def get_recal_input(bai=False):
#     # case 1: no duplicate removal
#     f = "mapped/{sample}-{unit}.sorted.bam"
#     if config["processing"]["remove-duplicates"]:
#         # case 2: remove duplicates
#         f = "dedup/{sample}-{unit}.bam"
#     if bai:
#         if config["processing"].get("restrict-regions"):
#             # case 3: need an index because random access is required
#             f += ".bai"
#             return f
#         else:
#             # case 4: no index needed
#             return []
#     else:
#         return f


########################################################################################################################
## Load genome and build index for alignment and gatk
path_genome = str(path_data + "genome/" + config["ref"]["name"] + "/" + config["ref"]["genome"].split("/")[-1])
path_dict = path_genome.replace("fa", "dict").replace("fasta", "dict")
path_dict_orginal = config["ref"]["genome"].replace("fa", "dict").replace("fasta", "dict")
path_genome_list=[ path_genome+".fai",path_genome+".bwt",path_dict,path_dict+"_chrom"]
rule LoadGenome:
    input:
         config["ref"]["genome"]
    output:
          path_genome
    shell:
         """
         ln -sr {input} {output}
         touch -h {output}
         """
rule GenomeIndexSamtools:
    input:
         path_genome
    output:
          path_genome + ".fai",  # samtools faidx
    log:
       path_log + "genomeindex/genomeindex_samtools.log"
    run:
        if os.path.exists(config["ref"]["genome"] + ".fai"):
            shell("ln -sr {orign_ref}.fai {path_genome}.fai".format(orign_ref=config["ref"]["genome"],
                                                                    path_genome=path_genome))
        else:
            shell("{path_samtools}samtools faidx {input}  2>{log} 1>{log} ")

rule GenomeIndexBwa:
    input:
         path_genome
    output:
          path_genome + ".amb",  # for bwa
          path_genome + ".bwt",  # for bwa
          path_genome + ".pac",  # for bwa
          path_genome + ".ann",  # for bwa
          path_genome + ".sa",  # for bwa
    log:
       path_log + "genomeindex/genomeindex_bwa.log"
    run:
        bwt_suffix = ["amb", "ann", "bwt", "pac", "sa"]
        bwt_stat = []
        for i in bwt_suffix:
            bwt_stat.append(os.path.exists(config["ref"]["genome"] + "." + i))
        if (False in bwt_stat):
            shell("{path_bwa}bwa index {input}  2>{log} 1>{log}")
        else:
            for i in bwt_suffix:
                shell("ln -sr " + config["ref"]["genome"] + "." + i + " " + path_genome + "." + i)
                shell("touch -h {path_genome}.{i}")
rule GenomeIndexPicard:
    input:
         path_genome
    output:
          path_dict  # for gatk
    log:
       path_log + "genomeindex/genomeindex_picard.log"
    run:
        if os.path.exists(path_dict_orginal):
            shell("ln -sr {path_dict_orginal} {path_dict}")
        else:
            shell("{path_picard}picard CreateSequenceDictionary R={path_genome} O={path_dict}  2>{log} 1>{log} ")
