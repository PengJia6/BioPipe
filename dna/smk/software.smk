########################################################################################################################
### define software path
default = "/home/pengjia/miniconda3/envs/ngs/bin/"
default27 = "/home/pengjia/miniconda3/envs/ngs27/bin/"
default = default if default[-1] == "/" else default + "/"
default27 = default27 if default27[-1] == "/" else default27 + "/"
mainEnv = "/home/pengjia/miniconda3/envs/PQsnake/bin/"
# qc.smk

path_fastqc = default
path_fastp = default
path_multiqc = default
path_pigz = default
path_samtools = default
path_bwa = default
path_picard = default
path_biobambam = default
path_bcftools = default
path_bedtools = default
path_gatk3 = default  # conda install and gatk3-register; the GenomeAnalysisTK.jar is jar package with version 3.8
path_gatk = default
path_freebayes = default
path_varscan = default
path_bgzip = default
path_vcflib = default
path_tabix = default

path_cnvnator = default
path_breakdancer = default
path_delly = default


path_breakseq2=default27
path_lumpy_sv = default27
path_smoove = default27
path_metasv = default27
path_spades = default27
path_age = "/home/pengjia/mysoftware/sv/AGE/age_v0.4_for_metasv/src/"
