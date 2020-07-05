path_raw_vcf = []
if "HC" in config["pipe"]["snpindel"]:
    # single sample calling
    for i in bam_sample_list:
        path_raw_vcf.append(
            path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.raw.vcf.gz".format(bam_sample=i))
if "HC_Hard_filter" in config["pipe"]["snpindel"]:
    for i in bam_sample_list:
        path_raw_vcf.append(
            path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.SNV.passh.vcf.gz".format(bam_sample=i))
        path_raw_vcf.append(
            path_data + "germlineVar/HC/perSample/{bam_sample}/{bam_sample}.HC.INDEL.passh.vcf.gz".format(bam_sample=i))
if "HC_Joint" in config["pipe"]["snpindel"]:
    # joint calling
    path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.raw.vcf.gz")

if "HC_Jonit_Hard_filter" in config["pipe"]["snpindel"]:
     path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.SNV.passh.vcf.gz")
     path_raw_vcf.append(path_data + "germlineVar/HC/jointCall/jointCall/" + config["project"]["name"] + ".HC.INDEL.passh.vcf.gz")

if "freebayes" in config["pipe"]["snpindel"]:
    for i in bam_sample_list:
        path_raw_vcf.append(
            path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.raw.vcf.gz.tbi".format(bam_sample=i))
        path_raw_vcf.append(
            path_data + "germlineVar/FB/perSample/{bam_sample}/{bam_sample}.FB.pass.vcf.gz.tbi".format(bam_sample=i))
if "samtools" in config["pipe"]["snpindel"]:
    for i in bam_sample_list:
        path_raw_vcf.append(
            path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.raw.vcf.gz.tbi".format(bam_sample=i))
        path_raw_vcf.append(
            path_data + "germlineVar/Samtools/perSample/{bam_sample}/{bam_sample}.Samtools.pass.vcf.gz.tbi".format(bam_sample=i))

