import cyvcf2

vcf_path = "/re_gecip/shared_allGeCIPs/vchundru/DNMs/v2_aggv3/aggv3_dnms_plusDecode_rmSegDupsLCRs.vcf.gz"

vcf = cyvcf2.VCF(vcf_path)

for variant in vcf:
    chrom = variant.CHROM
    pos = variant.POS
    ref = variant.REF
    alt = ",".join(variant.ALT)
    qual = variant.QUAL
    filter_pass = variant.FILTER
    info = variant.INFO

    print(f"{chrom}:{pos} {ref}>{alt} QUAL={qual} FILTER={filter_pass}")
    print("  INFO keys:", list(info.keys()))

    for sample in vcf.samples:
        gt = variant.genotypes[vcf.samples.index(sample)]
        # gt is [allele1, allele2, phased]
        print(f"  {sample} GT={gt[0]}/{gt[1]} phased={gt[2]}")

    # stop after the first variant if you only want a preview
    break