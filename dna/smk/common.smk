localrules: LoadGenome, LoadGenomeFile, GenomeIndexSamtools, GenomeIndexPicard
from snakemake.io import expand
import pandas as pd
import os
from snakemake.shell import shell

# report: "../reports/workflow.rst"
configfile: "conf/config.yaml"
path_data = os.path.abspath(config["path_data"]).rstrip("/") + "/"
path_log = os.path.abspath(config["path_log"]).rstrip("/") + "/"
path_bm = os.path.abspath(config["path_bm"]).rstrip("/") + "/"

########################################################################################################################
## Load genome and build index for alignment and gatk


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
# path_HQbamList = expand([path_data + "HQbam/{u.case}_{u.sample}.bam.bai"], u=caseinfo.itertuples())
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


#
# path_raw_vcf = []
# if "HC" in config["pipe"]["snpindel"]:
#     # single sample calling
#     for i in bam_sample_list:
#         path_raw_vcf.append(
#             path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.raw.vcf.gz".format(bam_sample=i))
# if "HC_Hard_filter" in config["pipe"]["snpindel"]:
#     for i in bam_sample_list:
#         path_raw_vcf.append(
#             path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.SNV.passh.vcf.gz".format(bam_sample=i))
#         path_raw_vcf.append(
#             path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.INDEL.passh.vcf.gz".format(bam_sample=i))
# if "HC_Joint" in config["pipe"]["snpindel"]:
#     # joint calling
#     path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.raw.vcf.gz")
#
# if "HC_Jonit_Hard_filter" in config["pipe"]["snpindel"]:
#      path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.SNV.passh.vcf.gz")
#      path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.INDEL.passh.vcf.gz")
#
# if "freebayes" in config["pipe"]["snpindel"]:
#     for i in bam_sample_list:
#         path_raw_vcf.append(
#             path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf.gz.tbi".format(bam_sample=i))
#         path_raw_vcf.append(
#             path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf.gz.tbi".format(bam_sample=i))
# if "samtools" in config["pipe"]["snpindel"]:
#     for i in bam_sample_list:
#         path_raw_vcf.append(
#             path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.raw.vcf.gz.tbi".format(bam_sample=i))
#         path_raw_vcf.append(
#             path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.pass.vcf.gz.tbi".format(bam_sample=i))
#

##### Wildcard constraints #####
wildcard_constraints:
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
                    bam_sample="|".join(bam_sample_list),
                    vartype="|".join(["SNV", "INDEL"])

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
        platform=caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), "PL"],
        LB=caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), "LB"])


# for bam_merge input
def get_sample_bam(wildcards):
    units = list(set(caseinfo.loc[(wildcards.case, wildcards.sample), "unit"]))
    aligner = wildcards.aligner
    case = wildcards.case
    sample = wildcards.sample
    qcpipe = wildcards.qcpipe
    return [path_data + "align/" + case + "/" + sample + "/" + case + "_" + sample + "_" + str(
        i) + "_" + qcpipe + "_" + aligner + "_sorted.bam" for i in units]
