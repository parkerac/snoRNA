#!/usr/bin/env python3
"""
Summarize snoRNA/scaRNA variants from a TSV/TSV.gz file and generate two summaries:
1) per-gene counts of unique variants and unique participants
2) per-variant output with small-RNA gene annotation

Usage:
  python 2_summarize_snoRNA_tsv.py \
    --tsv variants.tsv.gz \
    --gtf gencode.v49.annotation.gtf.gz \
    --out-gene-summary snoRNA_gene_summary.tsv \
    --out-variant-summary snoRNA_variant_summary.tsv

This script uses the standard library only.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

from snoRNA_utils import (
    DEFAULT_RNA_CLASSES,
    find_overlapping_snoRNAs,
    open_maybe_gzip,
    parse_feature_types,
    parse_snoRNA_regions,
    require_columns,
)


def summarize_tsv(tsv_path, gtf_path, out_gene_summary, out_variant_summary, feature_types=None):
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    feature_types = parse_feature_types(feature_types)

    print("Parsing GTF for small-RNA regions...", file=sys.stderr)
    regions_by_chrom = parse_snoRNA_regions(gtf_path, feature_types=feature_types)

    variant_counts = defaultdict(set)
    participant_counts = defaultdict(set)
    variant_rows = []

    with open_maybe_gzip(tsv_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {'ParticipantId', 'VariantId', 'chr', 'start', 'end'}, 'TSV')

        total_rows = 0
        filtered_rows = 0
        for row in reader:
            total_rows += 1
            try:
                chrom = row['chr']
                start = int(row['start'])
                end = int(row['end'])
            except Exception:
                continue
            if end < start:
                start, end = end, start
            rna_genes = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)
            if not rna_genes:
                continue

            filtered_rows += 1
            participant_id = row['ParticipantId']
            variant_id = row['VariantId']
            for rna_class, gene_name, gene_id in rna_genes:
                gene_key = (rna_class, gene_name, gene_id)
                variant_counts[gene_key].add(variant_id)
                participant_counts[gene_key].add(participant_id)
                output_row = dict(row)
                output_row['rna_class'] = rna_class
                output_row['rna_gene'] = gene_name
                output_row['snoRNA_gene'] = gene_name
                output_row['gene_id'] = gene_id
                variant_rows.append(output_row)

    print(f"Total rows read: {total_rows}", file=sys.stderr)
    print(f"Filtered small-RNA rows: {filtered_rows}", file=sys.stderr)
    print(f"Unique small-RNA-overlapping variants: {len({(r['VariantId'], r['chr'], r['start'], r['end'], r['ref'], r['alt']) for r in variant_rows})}", file=sys.stderr)

    print(f"Writing gene summary to {out_gene_summary}...", file=sys.stderr)
    with open(out_gene_summary, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['rna_class', 'gene_name', 'gene_id', 'unique_variant_count', 'unique_participant_count'])
        gene_summary = []
        for gene_key, variants in variant_counts.items():
            if not variants:
                continue
            rna_class, gene_name, gene_id = gene_key
            gene_summary.append((rna_class, gene_name, gene_id, len(variants), len(participant_counts[gene_key])))
        gene_summary.sort(key=lambda x: (-x[3], x[0], x[1], x[2]))
        for rna_class, gene_name, gene_id, variant_count, participant_count in gene_summary:
            writer.writerow([rna_class, gene_name, gene_id, variant_count, participant_count])

    print(f"Writing variant-level summary to {out_variant_summary}...", file=sys.stderr)
    with open(out_variant_summary, 'w', newline='') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=list(variant_rows[0].keys()) if variant_rows else ['ParticipantId', 'VariantId', 'chr', 'start', 'end', 'ref', 'alt', 'overlap_segdup_lcr', 'study', 'paternal_age', 'maternal_age', 'rna_class', 'rna_gene', 'snoRNA_gene', 'gene_id'])
        writer.writeheader()
        for row in variant_rows:
            writer.writerow(row)

    print(f"Gene summary contains {len(gene_summary)} small-RNA genes", file=sys.stderr)
    print(f"Variant-level summary contains {len(variant_rows)} rows", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Summarize snoRNA/scaRNA variant counts from TSV/TSV.gz source')
    parser.add_argument('--tsv', required=True, help='Input TSV or TSV.gz file with variant calls')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify small-RNA regions')
    parser.add_argument('--feature-types', default=','.join(DEFAULT_RNA_CLASSES), help='Comma-separated gene_type values to summarize (default: snoRNA,scaRNA)')
    parser.add_argument('--out-gene-summary', required=True, help='Output TSV path for small-RNA gene summary')
    parser.add_argument('--out-variant-summary', required=True, help='Output TSV path for variant-level summary')
    args = parser.parse_args()

    summarize_tsv(args.tsv, args.gtf, args.out_gene_summary, args.out_variant_summary, args.feature_types)


if __name__ == '__main__':
    main()
