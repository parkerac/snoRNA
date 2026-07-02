#!/usr/bin/env python3
"""
Summarize a filtered snoRNA VCF by counting unique variants per snoRNA gene.

Usage:
  python summarize_snoRNA_vcf.py --vcf snoRNA_only.vcf --gtf gencode.v49.annotation.gtf.gz --out summary.tsv

This script requires `cyvcf2`:
  python3 -m pip install cyvcf2

Output: TSV with columns: gene_name, chrom, start, end, unique_variant_count
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from bisect import bisect_right


def parse_gtf_for_snoRNA_regions(gtf_path):
    """Parse GTF to extract merged regions for all snoRNA features (any type with gene_type='snoRNA').
    Returns dict: chrom -> sorted list of (start, end, gene_name).
    Merges overlapping regions from the same gene.
    Accepts plain or gzipped GTFs.
    """
    regions_by_gene = defaultdict(lambda: defaultdict(list))  # gene_name -> chrom -> [(start, end), ...]
    
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, src, typ, start, end, score, strand, frame, attrs = parts
            
            # Parse attributes
            attr_dict = {}
            for kv in attrs.split(';'):
                kv = kv.strip()
                if not kv or ' ' not in kv:
                    continue
                k, v = kv.split(' ', 1)
                v = v.strip().strip('"')
                attr_dict[k] = v
            
            gene_type = attr_dict.get('gene_type') or attr_dict.get('gene_biotype')
            if gene_type != 'snoRNA':
                continue
            
            gene_name = attr_dict.get('gene_name', 'UNKNOWN')
            s = int(start)
            e = int(end)
            regions_by_gene[gene_name][chrom].append((s, e))
    
    # Merge overlapping regions per gene per chrom
    result = {}  # chrom -> [(start, end, gene_name), ...]
    default_dict = defaultdict(list)
    
    for gene_name, chroms in regions_by_gene.items():
        for chrom, intervals in chroms.items():
            intervals.sort()
            merged = []
            cur_s, cur_e = intervals[0]
            for s, e in intervals[1:]:
                if s <= cur_e + 1:
                    cur_e = max(cur_e, e)
                else:
                    merged.append((cur_s, cur_e, gene_name))
                    cur_s, cur_e = s, e
            merged.append((cur_s, cur_e, gene_name))
            default_dict[chrom].extend(merged)
    
    # Convert to regular dict and sort each chrom
    for chrom in default_dict:
        result[chrom] = sorted(default_dict[chrom])
    
    return result


def find_overlapping_genes(chrom, pos, regions_by_chrom):
    """Find all snoRNA genes overlapping the given position using binary search.
    regions_by_chrom: dict chrom -> sorted list of (start, end, gene_name)
    """
    if chrom not in regions_by_chrom:
        return []
    
    region_list = regions_by_chrom[chrom]
    overlapping = []
    
    # Binary search to find candidate regions
    starts = [s for s, e, g in region_list]
    idx = bisect_right(starts, pos)
    
    # Check the region just before and at idx
    for i in range(max(0, idx - 1), min(len(region_list), idx + 1)):
        s, e, gene_name = region_list[i]
        if s <= pos <= e:
            overlapping.append(gene_name)
    
    return overlapping


def summarize_vcf(vcf_path, gtf_path, out_path):
    """Count unique variants per snoRNA gene and write TSV summary."""
    try:
        from cyvcf2 import VCF
    except Exception as e:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise
    
    # Parse GTF for snoRNA regions
    print("Parsing GTF for snoRNA regions...", file=sys.stderr)
    regions_by_chrom = parse_gtf_for_snoRNA_regions(gtf_path)
    
    # Build gene_info from regions
    gene_info = {}  # gene_name -> (chrom, start, end)
    seen_genes = set()
    for chrom, region_list in regions_by_chrom.items():
        for start, end, gene_name in region_list:
            if gene_name not in seen_genes:
                gene_info[gene_name] = (chrom, start, end)
                seen_genes.add(gene_name)
    
    print(f"Found {len(gene_info)} snoRNA genes in GTF", file=sys.stderr)
    
    # Count variants per gene
    print("Reading VCF and counting variants...", file=sys.stderr)
    variant_counts = defaultdict(set)  # gene_name -> set of variant IDs
    vcf = VCF(vcf_path)
    
    total_variants = 0
    matched_variants = 0
    unmatched_chroms = set()
    
    for rec in vcf:
        total_variants += 1
        chrom = rec.CHROM
        pos = rec.POS
        variant_id = f"{chrom}:{pos}_{rec.REF}_{''.join(rec.ALT)}"
        
        overlapping = find_overlapping_genes(chrom, pos, regions_by_chrom)
        if not overlapping:
            # Variant doesn't overlap any snoRNA region
            if chrom not in regions_by_chrom:
                unmatched_chroms.add(chrom)
            continue
        
        matched_variants += 1
        for gene_name in overlapping:
            variant_counts[gene_name].add(variant_id)
    
    vcf.close()
    
    # Write output TSV (only genes with variants, sorted by count descending)
    print("Writing summary...", file=sys.stderr)
    with open(out_path, 'w') as out:
        out.write("gene_name\tchrom\tstart\tend\tunique_variant_count\n")
        
        # Filter to only genes with at least one variant and sort by count descending
        genes_with_variants = [
            (gene_name, len(variant_counts[gene_name]))
            for gene_name in gene_info.keys()
            if len(variant_counts[gene_name]) > 0
        ]
        genes_with_variants.sort(key=lambda x: x[1], reverse=True)
        
        for gene_name, count in genes_with_variants:
            chrom, start, end = gene_info[gene_name]
            out.write(f"{gene_name}\t{chrom}\t{start}\t{end}\t{count}\n")
    
    print(f"Summary written to {out_path}", file=sys.stderr)
    print(f"Total variants in VCF: {total_variants}", file=sys.stderr)
    print(f"Variants matched to snoRNA regions: {matched_variants}", file=sys.stderr)
    print(f"Total snoRNA genes: {len(gene_info)}", file=sys.stderr)
    print(f"Genes with variants: {sum(1 for c in variant_counts.values() if c)}", file=sys.stderr)
    if unmatched_chroms:
        print(f"WARNING: Chromosomes in VCF not found in GTF: {sorted(unmatched_chroms)}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description="Summarize snoRNA VCF by counting unique variants per gene")
    p.add_argument('--vcf', required=True, help='Filtered snoRNA VCF')
    p.add_argument('--gtf', required=True, help='GTF file (gz ok)')
    p.add_argument('--out', required=True, help='Output TSV path')
    args = p.parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(2)
    
    if not os.path.exists(args.gtf):
        print(f"GTF not found: {args.gtf}", file=sys.stderr)
        sys.exit(2)
    
    summarize_vcf(args.vcf, args.gtf, args.out)


if __name__ == '__main__':
    main()
