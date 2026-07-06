#!/usr/bin/env python3
"""
Summarize snoRNA/scaRNA variant carriers by NDD status with GEL denominators and deCODE separation.

This script reads a small-RNA variant TSV or TSV.gz file and a GEL phenotype TSV.
It categorizes each Aggv3 participant into:
  - undiagnosed NDD
  - diagnosed NDD
  - other

The output table is one row per RNA class / gene / gene_id and includes:
  - Aggv3 counts by group with GEL denominators in the same columns
  - deCODE variant carrier counts

Rows are sorted by RNA class, gene, and gene_id.

Usage:
  python 5_summarize_snoRNA_ndd_status.py \
    --variants variants.tsv.gz \
    --phenotype phenotype.tsv \
    --gtf gencode.v49.annotation.gtf.gz \
    --out snoRNA_ndd_status_summary.tsv

If the variant file already contains RNA annotation columns
(`rna_class`, `rna_gene`/`snoRNA_gene`, and optionally `gene_id`), the script
will use those values. Otherwise it uses the provided GTF to assign RNA genes
from variant intervals.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

from snoRNA_utils import (
    DEFAULT_RNA_CLASSES,
    aggv3_status_key,
    ensure_ndd_denominators,
    find_overlapping_snoRNAs,
    get_rna_annotation,
    load_ndd_phenotypes,
    open_maybe_gzip,
    parse_feature_types,
    parse_snoRNA_regions,
    require_columns,
)


def summarize_snoRNA_ndd_status(variants_path, phenotype_path, gtf_path, out_path, feature_types=None):
    if not os.path.exists(variants_path):
        raise FileNotFoundError(f'Variants file not found: {variants_path}')
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f'GTF file not found: {gtf_path}')

    print('Loading phenotype data...', file=sys.stderr)
    phenotypes, phenotype_denominators = load_ndd_phenotypes(phenotype_path)
    ensure_ndd_denominators(phenotype_denominators)

    feature_types = parse_feature_types(feature_types)

    print('Parsing GTF for small-RNA gene intervals...', file=sys.stderr)
    regions_by_chrom = parse_snoRNA_regions(gtf_path, feature_types=feature_types)

    gene_stats = defaultdict(lambda: {
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
            study = str(row['study']).strip()

            matches = []
            annotation = get_rna_annotation(row)
            if annotation:
                matches.append(annotation)
            else:
                try:
                    chrom = row['chr']
                    start = int(row['start'])
                    end = int(row['end'])
                except Exception:
                    continue
                matches = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)

            for rna_class, gene_name, gene_id in matches:
                summary = gene_stats[(rna_class, gene_name, gene_id)]
                if study.lower() == 'decode':
                    summary['deCODE_participants'].add(pid)
                    continue

                group = phenotypes.get(pid, 'other')
                summary[aggv3_status_key(group)].add(pid)

    print(f'Writing summary to {out_path}...', file=sys.stderr)
    with open(out_path, 'w', newline='') as out:
        fieldnames = [
            'rna_class',
            'snoRNA_gene',
            'rna_gene',
            'gene_id',
            'aggv3_undiagnosed_ndd',
            'aggv3_diagnosed_ndd',
            'aggv3_other',
            'deCODE_variant_carriers',
        ]
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for (rna_class, gene_name, gene_id), stats in sorted(
            gene_stats.items(),
            key=lambda x: (x[0][0], x[0][1], x[0][2]),
        ):
            writer.writerow({
                'rna_class': rna_class,
                'snoRNA_gene': gene_name,
                'rna_gene': gene_name,
                'gene_id': gene_id,
                'aggv3_undiagnosed_ndd': f"{len(stats['aggv3_undiagnosed'])}/{phenotype_denominators['undiagnosed NDD']}",
                'aggv3_diagnosed_ndd': f"{len(stats['aggv3_diagnosed'])}/{phenotype_denominators['diagnosed NDD']}",
                'aggv3_other': f"{len(stats['aggv3_other'])}/{phenotype_denominators['other']}",
                'deCODE_variant_carriers': len(stats['deCODE_participants']),
            })
    print('Done.', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Summarize snoRNA/scaRNA NDD status with GEL denominators and deCODE separation')
    parser.add_argument('--variants', required=True, help='Input variant TSV or TSV.gz file with ParticipantId, VariantId, chr, start, end, study')
    parser.add_argument('--phenotype', required=True, help='GEL phenotype TSV file with platekey, ndd, and case_solved')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify small-RNA regions when variant annotation is missing')
    parser.add_argument('--feature-types', default=','.join(DEFAULT_RNA_CLASSES), help='Comma-separated gene_type values to summarize (default: snoRNA,scaRNA)')
    parser.add_argument('--out', required=True, help='Output TSV path')
    args = parser.parse_args()
    summarize_snoRNA_ndd_status(args.variants, args.phenotype, args.gtf, args.out, args.feature_types)


if __name__ == '__main__':
    main()
