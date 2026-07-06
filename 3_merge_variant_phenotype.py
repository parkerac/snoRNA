#!/usr/bin/env python3
"""
Merge small-RNA variant summary with phenotype data by participant ID.

Input:
  - variant summary TSV with ParticipantId, VariantId, and rna_gene/snoRNA_gene
  - phenotype TSV with platekey and hpo_terms

Output:
  - TSV with columns: rna_class, rna_gene, snoRNA, VariantId, ParticipantId, hpo_terms
  - sorted by RNA class, RNA gene, then VariantId

Usage:
  python 3_merge_variant_phenotype.py \
    --variant-summary snoRNA_variant_summary.tsv \
    --phenotype phenotype_data.tsv \
    --out snoRNA_variant_phenotype.tsv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

from snoRNA_utils import get_rna_class, get_rna_gene, require_columns


def load_phenotypes(phenotype_path):
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")
    data = defaultdict(set)
    with open(phenotype_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {'platekey', 'hpo_terms'}, 'Phenotype file')
        for row in reader:
            key = row['platekey']
            terms = row.get('hpo_terms', '').strip()
            if terms:
                data[key].add(terms)
            else:
                data[key].add('')
    # collapse duplicate entries
    return {k: ';'.join(sorted({t for t in vals if t})) if any(vals) else '' for k, vals in data.items()}


def merge_variant_phenotype(variant_path, phenotype_path, out_path):
    if not os.path.exists(variant_path):
        raise FileNotFoundError(f"Variant summary file not found: {variant_path}")
    phenotypes = load_phenotypes(phenotype_path)

    rows = []
    with open(variant_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {'ParticipantId', 'VariantId'}, 'Variant summary')
        if not {'rna_gene', 'snoRNA_gene', 'snoRNA'} & set(reader.fieldnames or []):
            raise ValueError('Variant summary must contain rna_gene, snoRNA_gene, or snoRNA column')
        for row in reader:
            participant_id = row['ParticipantId']
            variant_id = row['VariantId']
            rna_class = get_rna_class(row)
            rna_gene = get_rna_gene(row)
            hpo_terms = phenotypes.get(participant_id, '')
            rows.append({
                'rna_class': rna_class,
                'rna_gene': rna_gene,
                'snoRNA': rna_gene,
                'VariantId': variant_id,
                'ParticipantId': participant_id,
                'hpo_terms': hpo_terms,
            })

    rows.sort(key=lambda x: (x['rna_class'] or '', x['rna_gene'] or '', x['VariantId'] or ''))

    with open(out_path, 'w', newline='') as out:
        fieldnames = ['rna_class', 'rna_gene', 'snoRNA', 'VariantId', 'ParticipantId', 'hpo_terms']
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Wrote merged output to {out_path}", file=sys.stderr)
    print(f"Rows written: {len(rows)}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Merge small-RNA variant summary with phenotype data')
    parser.add_argument('--variant-summary', required=True, help='TSV variant summary file with ParticipantId, VariantId, and rna_gene/snoRNA_gene/snoRNA')
    parser.add_argument('--phenotype', required=True, help='TSV file with platekey and hpo_terms')
    parser.add_argument('--out', required=True, help='Output TSV path')
    args = parser.parse_args()
    merge_variant_phenotype(args.variant_summary, args.phenotype, args.out)


if __name__ == '__main__':
    main()
