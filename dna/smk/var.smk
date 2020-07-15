include: "var/var_gatk_HC.smk"
include: "var/var_freebayes.smk"
include: "var/var_bcftools.smk"
include: "var/var_varscan.smk"
include: "var/var_delly.smk"
include: "var/var_smoove.smk"
include: "var/var_metasv.smk"
rule tabix:
    input:
         "{prefix}.vcf.gz"
    output:
          "{prefix}.vcf.gz.tbi"
    run:
        shell("{path_tabix}tabix -f {input} ")
