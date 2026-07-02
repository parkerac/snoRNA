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


def parse_gtf_for_snoRNA_genes(gtf_path):
    """Return list of (chrom, start, end, gene_name) for all snoRNA genes.
    Accepts plain or gzipped GTFs.
    """
    genes = []
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, src, typ, start, end, score, strand, frame, attrs = parts
            
            # Only process gene features
            if typ != 'gene':
                continue
            
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
            genes.append((chrom, s, e, gene_name))
    
    return genes


def find_overlapping_genes(chrom, pos, genes_by_chrom):
    """Find all snoRNA genes overlapping the given position.
    genes_by_chrom: dict chrom -> sorted list of (start, end, gene_name)
    """
    if chrom not in genes_by_chrom:
        return []
    
    gene_list = genes_by_chrom[chrom]
    overlapping = []
    for start, end, gene_name in gene_list:
        if start <= pos <= end:
            overlapping.append(gene_name)
    
    return overlapping


def summarize_vcf(vcf_path, gtf_path, out_path):
    """Count unique variants per snoRNA gene and write TSV summary."""
    try:
        from cyvcf2 import VCF
    except Exception as e:
        print("Missing dependency: cyvcf2 is required. Install with: python3 -m pip install cyvcf2", file=sys.stderr)
        raise
    
    # Parse GTF and organize by chrom
    print("Parsing GTF for snoRNA genes...", file=sys.stderr)
    genes = parse_gtf_for_snoRNA_genes(gtf_path)
    genes_by_chrom = defaultdict(list)
    gene_info = {}  # gene_name -> (chrom, start, end)
    
    for chrom, start, end, gene_name in genes:
        genes_by_chrom[chrom].append((start, end, gene_name))
        if gene_name not in gene_info:
            gene_info[gene_name] = (chrom, start, end)
    
    # Sort each chrom's genes by start
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort()
    
    # Count variants per gene
    print("Reading VCF and counting variants...", file=sys.stderr)
    variant_counts = defaultdict(set)  # gene_name -> set of variant IDs
    vcf = VCF(vcf_path)
    
    for rec in vcf:
        chrom = rec.CHROM
        pos = rec.POS
        variant_id = f"{chrom}:{pos}_{rec.REF}_{''.join(rec.ALT)}"
        
        overlapping = find_overlapping_genes(chrom, pos, genes_by_chrom)
        if not overlapping:
            # Variant doesn't overlap any snoRNA (shouldn't happen if VCF is filtered)
            continue
        
        for gene_name in overlapping:
            variant_counts[gene_name].add(variant_id)
    
    vcf.close()
    
    # Write output TSV
    print("Writing summary...", file=sys.stderr)
    with open(out_path, 'w') as out:
        out.write("gene_name\tchrom\tstart\tend\tunique_variant_count\n")
        
        for gene_name in sorted(gene_info.keys()):
            chrom, start, end = gene_info[gene_name]
            count = len(variant_counts[gene_name])
            out.write(f"{gene_name}\t{chrom}\t{start}\t{end}\t{count}\n")
    
    print(f"Summary written to {out_path}", file=sys.stderr)
    print(f"Total snoRNA genes: {len(gene_info)}", file=sys.stderr)
    print(f"Genes with variants: {sum(1 for c in variant_counts.values() if c)}", file=sys.stderr)


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
