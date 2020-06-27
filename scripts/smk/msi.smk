rule msisensor: 
    input:
      unpack(get_paired_bam)
    output:
      "../data/mutation/msi/msisnesor/{case}/{case}",
    log:
      "../logs/mutation/msi/msisensor/{case}.logs"
    benchmark:
      "../benchmark/mutation/msi/msisensor/{case}.tsv"
    params:
      ref="../data/common/GRCh38_full_analysis_set_plus_decoy_hla.fa.list",
      extra="",
      path=config["mainEnv"]
    threads: config["threads"]["msisensor"]
    wrapper:
      config["wrapper"]+"msisensor"

rule msisensor_pro: 
    input:
      unpack(get_paired_bam)
    output:
      N="../data/mutation/msi/msisnesor-pro-2/{case}/{case}-N",
      T="../data/mutation/msi/msisnesor-pro-2/{case}/{case}-T",
    log:
      "../logs/mutation/msi/msisensor-pro/{case}.logs"
    benchmark:
      "../benchmark/mutation/msi/msisensor-pro/{case}.tsv"
    params:
      ref="../data/common/GRCh38_full_analysis_set_plus_decoy_hla.fa.list",
      extra=" -i 0.2 ",
      path=config["mainEnv"]
    threads: config["threads"]["msisensor"]
    wrapper:
      config["wrapper"]+"msisensor-pro"

