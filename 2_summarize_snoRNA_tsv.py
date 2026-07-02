#!/usr/bin/env python3
"""
Summarize snoRNA variants from a TSV.gz file and generate two summaries:
1) per-snoRNA counts of unique variants and unique participants
2) per-variant output with snoRNA gene annotation

Usage:
  python 2_summarize_snoRNA_tsv.py \
    --tsv variants.tsv.gz \
    --gtf gencode.v49.annotation.gtf.gz \
    --out-gene-summary snoRNA_gene_summary.tsv \
    --out-variant-summary snoRNA_variant_summary.tsv

This script requires `cyvcf2` not used here; it uses `gzip` and the standard library.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from bisect import bisect_right


def parse_gtf_for_snoRNA_regions(gtf_path):
    """Parse GTF and return snoRNA regions merged by gene.

    Returns dict chrom -> sorted list of (start, end, gene_name).
    """
    regions_by_gene = defaultdict(lambda: defaultdict(list))
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, src, typ, start, end, score, strand, frame, attrs = parts
            attr_dict = {}
            for kv in attrs.split(';'):
                kv = kv.strip()
                if not kv or ' ' not in kv:
                    continue
                key, value = kv.split(' ', 1)
                attr_dict[key] = value.strip().strip('"')
            gene_type = attr_dict.get('gene_type') or attr_dict.get('gene_biotype')
            if gene_type != 'snoRNA':
                continue
            gene_name = attr_dict.get('gene_name', 'UNKNOWN')
            s = int(start)
            e = int(end)
            regions_by_gene[gene_name][chrom].append((s, e))

    regions_by_chrom = defaultdict(list)
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
            regions_by_chrom[chrom].extend(merged)

    for chrom in regions_by_chrom:
        regions_by_chrom[chrom].sort()
    return regions_by_chrom


def find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom):
    """Return snoRNA gene names overlapping the variant interval."""
    if chrom not in regions_by_chrom:
        return []
    regions = regions_by_chrom[chrom]
    starts = [s for s, e, g in regions]
    idx = bisect_right(starts, end)
    overlapping = []
    # Check previous regions and forward regions while they might overlap
    i = max(0, idx - 1)
    while i < len(regions) and regions[i][0] <= end:
        s, e, gene_name = regions[i]
        if s <= end and start <= e:
            overlapping.append(gene_name)
        i += 1
    return overlapping


def summarize_tsv(tsv_path, gtf_path, out_gene_summary, out_variant_summary):
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    print("Parsing GTF for snoRNA regions...", file=sys.stderr)
    regions_by_chrom = parse_gtf_for_snoRNA_regions(gtf_path)

    variant_counts = defaultdict(set)
    participant_counts = defaultdict(set)
    variant_rows = []

    with gzip.open(tsv_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'ParticipantId' not in reader.fieldnames:
            raise ValueError('TSV missing ParticipantId column')
        if 'VariantId' not in reader.fieldnames:
            raise ValueError('TSV missing VariantId column')
        if 'chr' not in reader.fieldnames or 'start' not in reader.fieldnames or 'end' not in reader.fieldnames:
            raise ValueError('TSV must include chr, start, and end columns')

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
            snoRNA_genes = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)
            if not snoRNA_genes:
                continue

            filtered_rows += 1
            participant_id = row['ParticipantId']
            variant_id = row['VariantId']
            for gene_name in snoRNA_genes:
                variant_counts[gene_name].add(variant_id)
                participant_counts[gene_name].add(participant_id)
                output_row = dict(row)
                output_row['snoRNA_gene'] = gene_name
                variant_rows.append(output_row)

    print(f"Total rows read: {total_rows}", file=sys.stderr)
    print(f"Filtered snoRNA rows: {filtered_rows}", file=sys.stderr)
    print(f"Unique snoRNA-overlapping variants: {len({(r['VariantId'], r['chr'], r['start'], r['end'], r['ref'], r['alt']) for r in variant_rows})}", file=sys.stderr)

    print(f"Writing gene summary to {out_gene_summary}...", file=sys.stderr)
    with open(out_gene_summary, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['gene_name', 'unique_variant_count', 'unique_participant_count'])
        gene_summary = [
            (gene_name, len(variant_counts[gene_name]), len(participant_counts[gene_name]))
            for gene_name in variant_counts.keys()
            if len(variant_counts[gene_name]) > 0
        ]
        gene_summary.sort(key=lambda x: x[1], reverse=True)
        for gene_name, variant_count, participant_count in gene_summary:
            writer.writerow([gene_name, variant_count, participant_count])

    print(f"Writing variant-level summary to {out_variant_summary}...", file=sys.stderr)
    with open(out_variant_summary, 'w', newline='') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=list(variant_rows[0].keys()) if variant_rows else ['ParticipantId', 'VariantId', 'chr', 'start', 'end', 'ref', 'alt', 'overlap_segdup_lcr', 'study', 'paternal_age', 'maternal_age', 'snoRNA_gene'])
        writer.writeheader()
        for row in variant_rows:
            writer.writerow(row)

    print(f"Gene summary contains {len(gene_summary)} snoRNAs", file=sys.stderr)
    print(f"Variant-level summary contains {len(variant_rows)} rows", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Summarize snoRNA variant counts from TSV.gz source')
    parser.add_argument('--tsv', required=True, help='Input TSV.gz file with variant calls')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify snoRNA regions')
    parser.add_argument('--out-gene-summary', required=True, help='Output TSV path for snoRNA gene summary')
    parser.add_argument('--out-variant-summary', required=True, help='Output TSV path for variant-level summary')
    args = parser.parse_args()

    summarize_tsv(args.tsv, args.gtf, args.out_gene_summary, args.out_variant_summary)


if __name__ == '__main__':
    main()
