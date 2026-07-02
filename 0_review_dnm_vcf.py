#!/usr/bin/env python3
"""
Filter a VCF to variants overlapping snoRNA features from a GTF.

Usage:
  python review_dnm_vcf.py --vcf /home/vscode/session_data/mounted-data-readonly/aggv3_dnms_plusDecode_rmSegDupsLCRs.vcf.gz --out snoRNA_only.vcf \
       [--gtf /path/to/gencode.v49.annotation.gtf.gz] [--feature snoRNA]

This script requires `cyvcf2` to read/write VCFs. Install with:
  python3 -m pip install cyvcf2

"""

import argparse
import gzip
import os
import sys
from bisect import bisect_right


def parse_gtf_for_feature(gtf_path, feature_type="snoRNA"):
    """Return dict chrom -> merged list of (start,end) intervals (1-based inclusive)
    for gene entries whose gene_type (attribute) equals feature_type.
    Accepts plain or gzipped GTFs.
    """
    regions = {}
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, src, typ, start, end, score, strand, frame, attrs = parts
            # we care about gene lines (or any line) that carry gene_type
            attr_dict = {}
            for kv in attrs.split(';'):
                kv = kv.strip()
                if not kv:
                    continue
                if ' ' in kv:
                    k, v = kv.split(' ', 1)
                    v = v.strip().strip('"')
                    attr_dict[k] = v
            gene_type = attr_dict.get('gene_type') or attr_dict.get('gene_biotype')
            if gene_type != feature_type:
                continue
            s = int(start)
            e = int(end)
            regions.setdefault(chrom, []).append((s, e))

    # merge intervals per chrom
    merged = {}
    for chrom, ivs in regions.items():
        ivs.sort()
        out = []
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                out.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        out.append((cur_s, cur_e))
        merged[chrom] = out
    return merged


def pos_overlaps_intervals(pos, ivs):
    """Return True if 1-based position `pos` overlaps any interval in sorted list `ivs`.
    `ivs` is list of (start,end) sorted by start.
    Uses binary search on starts.
    """
    if not ivs:
        return False
    starts = [s for s, e in ivs]
    i = bisect_right(starts, pos)
    # candidate interval is at i-1
    if i:
        s, e = ivs[i - 1]
        if s <= pos <= e:
            return True
    return False


def filter_vcf(vcf_in, vcf_out, gtf_path, feature="snoRNA"):
    try:
        from cyvcf2 import VCF, Writer
    except Exception as e:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise

    regions = parse_gtf_for_feature(gtf_path, feature_type=feature)
    vcf = VCF(vcf_in)
    w = Writer(vcf_out, vcf)
    count_in = 0
    count_out = 0
    for rec in vcf:
        count_in += 1
        chrom = rec.CHROM
        pos = rec.POS  # 1-based
        if pos_overlaps_intervals(pos, regions.get(chrom, [])):
            w.write_record(rec)
            count_out += 1
    w.close()
    vcf.close()
    print(f"Seen {count_in} variants; wrote {count_out} {feature} variants to {vcf_out}")


def main():
    p = argparse.ArgumentParser(description="Filter VCF for variants overlapping snoRNA features from GTF")
    p.add_argument('--vcf', required=True, help='Input VCF (can be .vcf or .vcf.gz)')
    p.add_argument('--out', required=True, help='Output VCF path')
    p.add_argument('--gtf', required=True, help='GTF path (gz ok)')
    p.add_argument('--feature', default='snoRNA', help='gene_type to filter for (default: snoRNA)')
    args = p.parse_args()

    if not os.path.exists(args.gtf):
        print(f"GTF not found at {args.gtf}. Provide --gtf to a valid GTF file.", file=sys.stderr)
        sys.exit(2)

    filter_vcf(args.vcf, args.out, args.gtf, feature=args.feature)


if __name__ == '__main__':
    main()

