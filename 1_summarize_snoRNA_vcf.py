#!/usr/bin/env python3
"""
Summarize a filtered small-RNA VCF by counting unique variants per gene.

Usage:
  python summarize_snoRNA_vcf.py --vcf small_rna_only.vcf --gtf gencode.v49.annotation.gtf.gz --out summary.tsv

This script requires `cyvcf2`:
  python3 -m pip install cyvcf2

Output: TSV with columns: rna_class, gene_name, gene_id, chrom, start, end, unique_variant_count
"""

import argparse
import os
import sys
from collections import defaultdict

from snoRNA_utils import DEFAULT_RNA_CLASSES, find_overlapping_snoRNAs, parse_feature_types, parse_snoRNA_regions


def summarize_vcf(vcf_path, gtf_path, out_path, feature_types=None):
    """Count unique variants per small-RNA gene and write TSV summary."""
    try:
        from cyvcf2 import VCF
    except Exception:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise
    
    feature_types = parse_feature_types(feature_types)

    print("Parsing GTF for small-RNA regions...", file=sys.stderr)
    regions_by_chrom = parse_snoRNA_regions(gtf_path, feature_types=feature_types)
    
    # Build gene_info from regions
    gene_info = {}  # (rna_class, gene_name, gene_id) -> (chrom, start, end)
    seen_genes = set()
    for chrom, region_list in regions_by_chrom.items():
        for start, end, rna_class, gene_name, gene_id in region_list:
            gene_key = (rna_class, gene_name, gene_id)
            if gene_key not in seen_genes:
                gene_info[gene_key] = (chrom, start, end)
                seen_genes.add(gene_key)
    
    print(f"Found {len(gene_info)} small-RNA genes in GTF", file=sys.stderr)
    
    # Count variants per gene
    print("Reading VCF and counting variants...", file=sys.stderr)
    variant_counts = defaultdict(set)  # (rna_class, gene_name, gene_id) -> set of variant IDs
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
            # Variant doesn't overlap any requested small-RNA region
            if chrom not in regions_by_chrom:
                unmatched_chroms.add(chrom)
            continue
        
        matched_variants += 1
        for rna_class, gene_name, gene_id in overlapping:
            variant_counts[(rna_class, gene_name, gene_id)].add(variant_id)
    
    vcf.close()
    
    # Write output TSV (only genes with variants, sorted by count descending)
    print("Writing summary...", file=sys.stderr)
    with open(out_path, 'w') as out:
        out.write("rna_class\tgene_name\tgene_id\tchrom\tstart\tend\tunique_variant_count\n")
        
        # Filter to only genes with at least one variant and sort by count descending
        genes_with_variants = [
            (gene_key, len(variant_counts[gene_key]))
            for gene_key in gene_info.keys()
            if len(variant_counts[gene_key]) > 0
        ]
        genes_with_variants.sort(key=lambda x: (-x[1], x[0]))
        
        for (rna_class, gene_name, gene_id), count in genes_with_variants:
            chrom, start, end = gene_info[(rna_class, gene_name, gene_id)]
            out.write(f"{rna_class}\t{gene_name}\t{gene_id}\t{chrom}\t{start}\t{end}\t{count}\n")
    
    print(f"Summary written to {out_path}", file=sys.stderr)
    print(f"Total variants in VCF: {total_variants}", file=sys.stderr)
    print(f"Variants matched to small-RNA regions: {matched_variants}", file=sys.stderr)
    print(f"Total small-RNA genes: {len(gene_info)}", file=sys.stderr)
    print(f"Genes with variants: {sum(1 for c in variant_counts.values() if c)}", file=sys.stderr)
    if unmatched_chroms:
        print(f"WARNING: Chromosomes in VCF not found in GTF: {sorted(unmatched_chroms)}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description="Summarize small-RNA VCF by counting unique variants per gene")
    p.add_argument('--vcf', required=True, help='Filtered small-RNA VCF')
    p.add_argument('--gtf', required=True, help='GTF file (gz ok)')
    p.add_argument('--out', required=True, help='Output TSV path')
    p.add_argument('--feature-types', default=','.join(DEFAULT_RNA_CLASSES), help='Comma-separated gene_type values to summarize (default: snoRNA,scaRNA)')
    args = p.parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"VCF not found: {args.vcf}", file=sys.stderr)
        sys.exit(2)
    
    if not os.path.exists(args.gtf):
        print(f"GTF not found: {args.gtf}", file=sys.stderr)
        sys.exit(2)
    
    summarize_vcf(args.vcf, args.gtf, args.out, feature_types=args.feature_types)


if __name__ == '__main__':
    main()
