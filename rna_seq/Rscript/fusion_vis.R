# Title     : TODO
# Objective : TODO
# Created by: pengjia
# Created on: 2021/1/9
library(chimeraviz)
library(chimeraviz)

#BiocManager::install("chimeraviz")
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19")

soapfuse833ke <- system.file(
  "extdata",
  "/mnt/project/ProjectSnake/star_fusion/star-fusion.fusion_predictions.tsv",
  package = "chimeraviz")
fusions <- import_starfusion(soapfuse833ke, "hg38")
fusions <- import_starfusion("/mnt/project/ProjectSnake/star_fusion/star-fusion.fusion_predictions.tsv",
                             "hg38",
                             limit = 10)
# Plot!
fusions
plot_circle(fusions)


fusion <- get_fusion_by_gene_name(fusions, "MIR205HG")
fusion

if (!exists("bamfile5267"))
  bamfile5267 <- system.file(
    "extdata",
    "/mnt/project/ProjectSnake/star_fusion/test/finspector.junction_reads.bam",
    #"/mnt/project/ProjectSnake/star_fusion/1447609-C_L1_fastp_STAR_Aligned.sortedByCoord.out.bam",
    package = "chimeraviz")
fusion <- add_fusion_reads_alignment(fusion,
                                     bamfile5267,)


plot_fusion_reads(fusion)

