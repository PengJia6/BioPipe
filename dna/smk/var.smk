include: "var/var_gatk_HC.smk"
include: "var/var_freebayes.smk"
include: "var/var_bcftools.smk"
rule tabix:
    input:
         "{prefix}.vcf.gz"
    output:
          "{prefix}.vcf.gz.tbi"
    run:
        shell("{path_tabix}tabix -f {input} ")
