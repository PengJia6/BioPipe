import pandas as pd
from snakemake.utils import validate
report: "../reports/workflow.rst"

###### Config file and sample sheets #####
configfile: "conf/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

caseinfo = pd.read_csv(config["caseinfo"]).set_index(["case","sample","unit"], drop=False)
#sampleinfo= pd.read_csv(config["caseinfo"]).set_index(["case","sample"], drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")

#units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
caseinfo.index = caseinfo.index.set_levels([i.astype(str) for i in caseinfo.index.levels])  # enforce str in index
#validate(units, schema="../schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_table(config["ref"]["genome"] + ".fai2",
                        header=None, usecols=[0], squeeze=True, dtype=str)

casedict={}
for index,info in caseinfo.iterrows():
    case=index[0]
    sample=index[1]
    if case not in casedict:
        casedict[case]={}
    casedict[case][sample]=""
num=0
sampleinfo=pd.DataFrame()
for case in casedict:
    for sample in casedict[case]:
        num+=1
        sampleinfo.loc[num,"case"]=case
        sampleinfo.loc[num,"sample"]=sample
        sampleinfo.loc[num,"id"]=case+sample
sampleinfo=sampleinfo.set_index(["case","sample"])
sampleinfo.index = sampleinfo.index.set_levels([i.astype(str) for i in sampleinfo.index.levels])  # enforce str in index
for case in casedict:
    for sample in casedict[case]:
        num+=1
        sampleinfo.loc[(case,sample),"case"]=case
        sampleinfo.loc[(case,sample),"sample"]=sample


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    qcpipe="fastp|trim|passqc",
    aligner="bwa|bowtie",
    rmDup="pcdRmDup|bioRmDup",
    realign="reAlign|noReAlign",
    R="1|2",
    case="|".join(caseinfo["case"]),
    sample="|".join(caseinfo["sample"]),
    unit="|".join(caseinfo["unit"]),
    LB="|".join(caseinfo["LB"]),
    contig="|".join(contigs)


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = caseinfo.loc[(wildcards.case, wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"R1": fastqs.fq1, "R2": fastqs.fq2}
    return {"R1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{case}_{sample}\tSM:{case}_{sample}\tPL:{platform}\tLB:{LB}'".format(
        sample=wildcards.sample,
        case=wildcards.case,
        platform=caseinfo.loc[(wildcards.case,wildcards.sample, wildcards.unit), "PL"],
        LB=caseinfo.loc[(wildcards.case,wildcards.sample, wildcards.unit), "LB"])

def get_sample_bam(wildcards):
    units=list(set(caseinfo.loc[(wildcards.case,wildcards.sample),"unit"]))
    aligner=wildcards.aligner
    case=wildcards.case
    sample=wildcards.sample
    qcpipe=wildcards.qcpipe
    return [ "../data/align/"+aligner+"/"+case+"/"+sample+"/"+case+"_"+sample+"_"+str(i)+"_"+qcpipe+"_"+aligner+"_sorted.bam" for i in units]

#    return units
    return [str("../data/align/"+wildcards.aligner+"/"+wildcards.case+"/"+wildcards.sample+"/"+wildcards.case+"_"+wildcards.sample+"_"+unit+"_"+wildcards.qcpipe+"_"+wildcards.aligner+"_sorted.bam---") for unit in units]
#return expand(["../data/align/"+wildcards.aligner+"/"+wildcards.case+"/"+wildcards.sample+"/"+wildcards.case+"_"+wildcards.sample+"_{unitD}_"+wildcards.qcpipe+"_"+wildcards.aligner+"_sorted.bam" ],unitD=["1","2"])
#return expand(["../data/align/"+wildcards.aligner+"/"+wildcards.case+"/"+wildcards.sample+"/"+wildcards.case+"_"+wildcards.sample+"_{unitD}_"+wildcards.qcpipe+"_"+wildcards.aligner+"_sorted.bam" ],unitD=caseinfo.loc[(wildcards.case,wildcards.sample)].unit)



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


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f
