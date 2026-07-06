#!/usr/bin/env python3
"""
Filter a VCF to variants overlapping snoRNA/scaRNA features from a GTF.

Usage:
  python review_dnm_vcf.py --vcf /home/vscode/session_data/mounted-data-readonly/aggv3_dnms_plusDecode_rmSegDupsLCRs.vcf.gz --out small_rna_only.vcf \
       [--gtf /path/to/gencode.v49.annotation.gtf.gz] [--feature-types snoRNA,scaRNA]

This script requires `cyvcf2` to read/write VCFs. Install with:
  python3 -m pip install cyvcf2

"""

import argparse
import os
import sys

from snoRNA_utils import DEFAULT_RNA_CLASSES, find_overlapping_snoRNAs, parse_feature_types, parse_snoRNA_regions


def filter_vcf(vcf_in, vcf_out, gtf_path, feature_types=None):
    try:
        from cyvcf2 import VCF, Writer
    except Exception:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise

    feature_types = parse_feature_types(feature_types)
    regions_by_chrom = parse_snoRNA_regions(gtf_path, feature_types=feature_types)
    vcf = VCF(vcf_in)
    w = Writer(vcf_out, vcf)
    count_in = 0
    count_out = 0
    for rec in vcf:
        count_in += 1
        chrom = rec.CHROM
        pos = rec.POS  # 1-based
        if find_overlapping_snoRNAs(chrom, pos, pos, regions_by_chrom):
            w.write_record(rec)
            count_out += 1
    w.close()
    vcf.close()
    print(f"Seen {count_in} variants; wrote {count_out} {','.join(feature_types)} variants to {vcf_out}")


def main():
    p = argparse.ArgumentParser(description="Filter VCF for variants overlapping snoRNA/scaRNA features from GTF")
    p.add_argument('--vcf', required=True, help='Input VCF (can be .vcf or .vcf.gz)')
    p.add_argument('--out', required=True, help='Output VCF path')
    p.add_argument('--gtf', required=True, help='GTF path (gz ok)')
    p.add_argument('--feature-types', default=','.join(DEFAULT_RNA_CLASSES), help='Comma-separated gene_type values to filter for (default: snoRNA,scaRNA)')
    p.add_argument('--feature', help='Single gene_type to filter for; overrides --feature-types')
    args = p.parse_args()

    if not os.path.exists(args.gtf):
        print(f"GTF not found at {args.gtf}. Provide --gtf to a valid GTF file.", file=sys.stderr)
        sys.exit(2)

    filter_vcf(args.vcf, args.out, args.gtf, feature_types=args.feature or args.feature_types)


if __name__ == '__main__':
    main()
