#!/usr/bin/env python3
"""
Create a variant-level NDD status summary for U3, SNORD118, and SNORA70.

This script reads a snoRNA variant TSV or TSV.gz file and a GEL phenotype TSV.
It generates one row per variant for the selected snoRNA genes and reports:
  - Aggv3 counts by NDD group with GEL denominators
  - deCODE variant carrier counts

Usage:
  python 6_summarize_snoRNA_ndd_status_variant_level.py \
    --variants variants.tsv.gz \
    --phenotype phenotype.tsv \
    --gtf gencode.v49.annotation.gtf.gz \
    --out snoRNA_variant_level_ndd_status.tsv

Only variants overlapping U3, SNORD118, or SNORA70 are included.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from bisect import bisect_right

TARGET_GENES = {'U3', 'SNORD118', 'SNORA70'}


def open_maybe_gzip(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def normalize_bool(value):
    if value is None:
        return None
    value = str(value).strip().lower()
    if value in {'true', '1', 'yes', 'y', 't'}:
        return True
    if value in {'false', '0', 'no', 'n', 'f'}:
        return False
    return None


def normalize_case_solved(value):
    if value is None:
        return 'unknown'
    value = str(value).strip().lower()
    if value in {'yes', 'y'}:
        return 'yes'
    if value in {'no', 'n'}:
        return 'no'
    if value in {'partially', 'partial'}:
        return 'partially'
    if value in {'unknown', ''}:
        return 'unknown'
    return value


def categorize_ndd_status(ndd_value, case_solved_value):
    ndd = normalize_bool(ndd_value)
    solved = normalize_case_solved(case_solved_value)
    if ndd is True:
        if solved == 'yes':
            return 'diagnosed NDD'
        return 'undiagnosed NDD'
    return 'other'


def parse_gtf_for_snoRNA_regions(gtf_path):
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
            gene_id = attr_dict.get('gene_id', 'UNKNOWN')
            s = int(start)
            e = int(end)
            regions_by_gene[(gene_name, gene_id)][chrom].append((s, e))

    regions_by_chrom = defaultdict(list)
    for (gene_name, gene_id), chroms in regions_by_gene.items():
        for chrom, intervals in chroms.items():
            intervals.sort()
            merged = []
            cur_s, cur_e = intervals[0]
            for s, e in intervals[1:]:
                if s <= cur_e + 1:
                    cur_e = max(cur_e, e)
                else:
                    merged.append((cur_s, cur_e, gene_name, gene_id))
                    cur_s, cur_e = s, e
            merged.append((cur_s, cur_e, gene_name, gene_id))
            regions_by_chrom[chrom].extend(merged)

    for chrom in regions_by_chrom:
        regions_by_chrom[chrom].sort()
    return regions_by_chrom


def find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom):
    if chrom not in regions_by_chrom:
        return []
    regions = regions_by_chrom[chrom]
    starts = [s for s, e, g, gid in regions]
    idx = bisect_right(starts, end)
    overlapping = []
    i = max(0, idx - 1)
    while i < len(regions) and regions[i][0] <= end:
        s, e, gene_name, gene_id = regions[i]
        if s <= end and start <= e:
            overlapping.append((gene_name, gene_id))
        i += 1
    return overlapping


def load_phenotypes(phenotype_path, platekey_col='platekey', ndd_col='ndd', case_solved_col='case_solved'):
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')
    phenotypes = {}
    totals = defaultdict(set)
    with open(phenotype_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if platekey_col not in reader.fieldnames:
            raise ValueError(f'Phenotype file must contain column {platekey_col}')
        if ndd_col not in reader.fieldnames:
            raise ValueError(f'Phenotype file must contain column {ndd_col}')
        if case_solved_col not in reader.fieldnames:
            raise ValueError(f'Phenotype file must contain column {case_solved_col}')
        for row in reader:
            pid = row[platekey_col]
            group = categorize_ndd_status(row.get(ndd_col), row.get(case_solved_col))
            phenotypes[pid] = group
            totals[group].add(pid)
    return phenotypes, {group: len(pids) for group, pids in totals.items()}


def summarize_variant_level(variants_path, phenotype_path, gtf_path, out_path):
    if not os.path.exists(variants_path):
        raise FileNotFoundError(f'Variants file not found: {variants_path}')
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f'GTF file not found: {gtf_path}')

    print('Loading phenotype data...', file=sys.stderr)
    phenotypes, phenotype_denominators = load_phenotypes(phenotype_path)
    phenotype_denominators.setdefault('undiagnosed NDD', 0)
    phenotype_denominators.setdefault('diagnosed NDD', 0)
    phenotype_denominators.setdefault('other', 0)

    print('Parsing GTF for snoRNA gene intervals...', file=sys.stderr)
    regions_by_chrom = parse_gtf_for_snoRNA_regions(gtf_path)

    row_stats = defaultdict(lambda: {
        'aggv3_undiagnosed': set(),
        'aggv3_diagnosed': set(),
        'aggv3_other': set(),
        'deCODE_participants': set(),
        'chrom': None,
        'start': None,
        'end': None,
    })

    with open_maybe_gzip(variants_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'ParticipantId', 'VariantId', 'chr', 'start', 'end', 'study'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Variants file must contain columns: {", ".join(sorted(missing))}')

        for row in reader:
            pid = row['ParticipantId']
            variant_id = row['VariantId']
            study = str(row['study']).strip()
            gene_name = row.get('snoRNA_gene') or row.get('snoRNA')
            gene_id = row.get('gene_id')

            matches = []
            if gene_name and gene_name in TARGET_GENES and gene_id:
                matches.append((gene_name, gene_id))
            else:
                try:
                    chrom = row['chr']
                    start = int(row['start'])
                    end = int(row['end'])
                except Exception:
                    continue
                if end < start:
                    start, end = end, start
                overlaps = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)
                matches = [(g, gid) for g, gid in overlaps if g in TARGET_GENES]

            if not matches:
                continue

            for gene_name, gene_id in matches:
                key = (gene_name, gene_id, variant_id)
                if row_stats[key]['chrom'] is None:
                    row_stats[key]['chrom'] = row.get('chr')
                    try:
                        row_stats[key]['start'] = int(row['start'])
                        row_stats[key]['end'] = int(row['end'])
                    except Exception:
                        row_stats[key]['start'] = None
                        row_stats[key]['end'] = None
                if study.lower() == 'decode':
                    row_stats[key]['deCODE_participants'].add(pid)
                else:
                    group = phenotypes.get(pid, 'other')
                    if group == 'diagnosed NDD':
                        row_stats[key]['aggv3_diagnosed'].add(pid)
                    elif group == 'undiagnosed NDD':
                        row_stats[key]['aggv3_undiagnosed'].add(pid)
                    else:
                        row_stats[key]['aggv3_other'].add(pid)

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
        for (gene_name, gene_id, variant_id), stats in row_stats.items():
            gene_totals[gene_name] += (
                len(stats['aggv3_undiagnosed']) + len(stats['aggv3_diagnosed']) + len(stats['aggv3_other'])
            )

        def sort_key(item):
            (gene_name, gene_id, variant_id), stats = item
            gene_total = gene_totals[gene_name]
            return (-gene_total, gene_name, variant_id)

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
    parser = argparse.ArgumentParser(description='Create a variant-level NDD status summary for U3, SNORD118, and SNORA70')
    parser.add_argument('--variants', required=True, help='Input variant TSV or TSV.gz file with ParticipantId, VariantId, chr, start, end, study')
    parser.add_argument('--phenotype', required=True, help='GEL phenotype TSV file with platekey, ndd, and case_solved')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify snoRNA regions when variant annotation is missing')
    parser.add_argument('--out', required=True, help='Output TSV path')
    args = parser.parse_args()
    summarize_variant_level(args.variants, args.phenotype, args.gtf, args.out)


if __name__ == '__main__':
    main()
