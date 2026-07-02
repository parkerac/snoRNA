#!/usr/bin/env python3
"""
Create a variant-level NDD status summary for all snoRNA variants.

This script reads a snoRNA variant TSV or TSV.gz file and a GEL phenotype TSV.
It generates one row per variant for all snoRNA genes and reports:
  - Aggv3 counts by NDD group with GEL denominators
  - deCODE variant carrier counts

Usage:
  python 6_summarize_snoRNA_ndd_status_variant_level.py \
    --variants variants.tsv.gz \
    --phenotype phenotype.tsv \
    --gtf gencode.v49.annotation.gtf.gz \
    --out snoRNA_variant_level_ndd_status.tsv

All snoRNA variants are included.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

from snoRNA_utils import (
    aggv3_status_key,
    ensure_ndd_denominators,
    find_overlapping_snoRNAs,
    load_ndd_phenotypes,
    natural_key,
    open_maybe_gzip,
    parse_snoRNA_regions,
    require_columns,
)


def summarize_variant_level(variants_path, phenotype_path, gtf_path, out_path):
    if not os.path.exists(variants_path):
        raise FileNotFoundError(f'Variants file not found: {variants_path}')
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f'GTF file not found: {gtf_path}')

    print('Loading phenotype data...', file=sys.stderr)
    phenotypes, phenotype_denominators = load_ndd_phenotypes(phenotype_path)
    ensure_ndd_denominators(phenotype_denominators)

    print('Parsing GTF for snoRNA gene intervals...', file=sys.stderr)
    regions_by_chrom = parse_snoRNA_regions(gtf_path)

    row_stats = defaultdict(lambda: {
        'aggv3_undiagnosed': set(),
        'aggv3_diagnosed': set(),
        'aggv3_other': set(),
        'deCODE_participants': set(),
    })

    with open_maybe_gzip(variants_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {'ParticipantId', 'VariantId', 'chr', 'start', 'end', 'study'}, 'Variants file')

        for row in reader:
            pid = row['ParticipantId']
            variant_id = row['VariantId']
            study = str(row['study']).strip()
            gene_name = row.get('snoRNA_gene') or row.get('snoRNA')
            gene_id = row.get('gene_id')

            if gene_name and gene_id:
                matches = [(gene_name, gene_id)]
            else:
                try:
                    chrom = row['chr']
                    start = int(row['start'])
                    end = int(row['end'])
                except Exception:
                    continue
                matches = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)

            if not matches:
                continue

            for gene_name, gene_id in matches:
                key = (gene_name, gene_id, variant_id)
                if study.lower() == 'decode':
                    row_stats[key]['deCODE_participants'].add(pid)
                    continue

                group = phenotypes.get(pid, 'other')
                row_stats[key][aggv3_status_key(group)].add(pid)

    print(f'Writing variant-level summary to {out_path}...', file=sys.stderr)
    with open(out_path, 'w', newline='') as out:
        fieldnames = [
            'snoRNA_gene',
            'gene_id',
            'VariantId',
            'aggv3_undiagnosed_ndd',
            'aggv3_diagnosed_ndd',
            'aggv3_other',
            'deCODE_variant_carriers',
        ]
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        gene_totals = defaultdict(int)
        for (gene_name, _gene_id, _variant_id), stats in row_stats.items():
            gene_totals[gene_name] += (
                len(stats['aggv3_undiagnosed']) + len(stats['aggv3_diagnosed']) + len(stats['aggv3_other'])
            )

        def sort_key(item):
            (gene_name, _gene_id, variant_id), _stats = item
            gene_total = gene_totals[gene_name]
            return (-gene_total, gene_name, natural_key(variant_id))

        for (gene_name, gene_id, variant_id), stats in sorted(row_stats.items(), key=sort_key):
            writer.writerow({
                'snoRNA_gene': gene_name,
                'gene_id': gene_id,
                'VariantId': variant_id,
                'aggv3_undiagnosed_ndd': f"{len(stats['aggv3_undiagnosed'])}/{phenotype_denominators['undiagnosed NDD']}",
                'aggv3_diagnosed_ndd': f"{len(stats['aggv3_diagnosed'])}/{phenotype_denominators['diagnosed NDD']}",
                'aggv3_other': f"{len(stats['aggv3_other'])}/{phenotype_denominators['other']}",
                'deCODE_variant_carriers': len(stats['deCODE_participants']),
            })
    print('Done.', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Create a variant-level NDD status summary for snoRNA variants')
    parser.add_argument('--variants', required=True, help='Input variant TSV or TSV.gz file with ParticipantId, VariantId, chr, start, end, study')
    parser.add_argument('--phenotype', required=True, help='GEL phenotype TSV file with platekey, ndd, and case_solved')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify snoRNA regions when variant annotation is missing')
    parser.add_argument('--out', required=True, help='Output TSV path')
    args = parser.parse_args()
    summarize_variant_level(args.variants, args.phenotype, args.gtf, args.out)


if __name__ == '__main__':
    main()
