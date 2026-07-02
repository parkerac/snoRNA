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
import os
import sys
from collections import defaultdict

from snoRNA_utils import find_overlapping_snoRNAs, parse_snoRNA_regions


def summarize_vcf(vcf_path, gtf_path, out_path):
    """Count unique variants per snoRNA gene and write TSV summary."""
    try:
        from cyvcf2 import VCF
    except Exception:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise
    
    # Parse GTF for snoRNA regions
    print("Parsing GTF for snoRNA regions...", file=sys.stderr)
    regions_by_chrom = parse_snoRNA_regions(gtf_path)
    
    # Build gene_info from regions
    gene_info = {}  # gene_name -> (chrom, start, end)
    seen_genes = set()
    for chrom, region_list in regions_by_chrom.items():
        for start, end, gene_name, _gene_id in region_list:
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
        
        overlapping = find_overlapping_snoRNAs(chrom, pos, pos, regions_by_chrom)
        if not overlapping:
            # Variant doesn't overlap any snoRNA region
            if chrom not in regions_by_chrom:
                unmatched_chroms.add(chrom)
            continue
        
        matched_variants += 1
        for gene_name, _gene_id in overlapping:
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
