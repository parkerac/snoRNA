#!/usr/bin/env python3
"""
Summarize snoRNA variant carriers by NDD status with GEL denominators and deCODE separation.

This script reads a snoRNA variant TSV or TSV.gz file and a GEL phenotype TSV.
It categorizes each Aggv3 participant into:
  - undiagnosed NDD
  - diagnosed NDD
  - other

The output table is one row per snoRNA gene / gene_id and includes:
  - Aggv3 counts by group with GEL denominators in the same columns
  - deCODE variant carrier counts

Usage:
  python 5_summarize_snoRNA_ndd_status.py \
    --variants variants.tsv.gz \
    --phenotype phenotype.tsv \
    --gtf gencode.v49.annotation.gtf.gz \
    --out snoRNA_ndd_status_summary.tsv

If the variant file already contains snoRNA gene annotation columns
(`snoRNA_gene` and optionally `gene_id`), the script will use those values.
Otherwise it uses the provided GTF to assign snoRNA genes from variant intervals.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict
from bisect import bisect_right


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


def summarize_snoRNA_ndd_status(variants_path, phenotype_path, gtf_path, out_path):
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

    gene_stats = defaultdict(lambda: {
        'aggv3_undiagnosed': set(),
        'aggv3_diagnosed': set(),
        'aggv3_other': set(),
        'deCODE_participants': set(),
        'deCODE_variants': set(),
        'total_aggv3_participants': set(),
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

            if not gene_name or not gene_id:
                try:
                    chrom = row['chr']
                    start = int(row['start'])
                    end = int(row['end'])
                except Exception:
                    continue
                if end < start:
                    start, end = end, start
                overlaps = find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom)
                if not overlaps:
                    continue
                for gene_name, gene_id in overlaps:
                    summary = gene_stats[(gene_name, gene_id)]
                    if study.lower() == 'decode':
                        summary['deCODE_participants'].add(pid)
                        summary['deCODE_variants'].add(variant_id)
                    else:
                        group = phenotypes.get(pid, 'other')
                        if group == 'diagnosed NDD':
                            summary['aggv3_diagnosed'].add(pid)
                        elif group == 'undiagnosed NDD':
                            summary['aggv3_undiagnosed'].add(pid)
                        else:
                            summary['aggv3_other'].add(pid)
                        summary['total_aggv3_participants'].add(pid)
            else:
                if not gene_id:
                    gene_id = row.get('gene_id', 'UNKNOWN')
                summary = gene_stats[(gene_name, gene_id)]
                if study.lower() == 'decode':
                    summary['deCODE_participants'].add(pid)
                    summary['deCODE_variants'].add(variant_id)
                else:
                    group = phenotypes.get(pid, 'other')
                    if group == 'diagnosed NDD':
                        summary['aggv3_diagnosed'].add(pid)
                    elif group == 'undiagnosed NDD':
                        summary['aggv3_undiagnosed'].add(pid)
                    else:
                        summary['aggv3_other'].add(pid)
                    summary['total_aggv3_participants'].add(pid)

    print(f'Writing summary to {out_path}...', file=sys.stderr)
    with open(out_path, 'w', newline='') as out:
        fieldnames = [
            'snoRNA_gene',
            'gene_id',
            'aggv3_undiagnosed_ndd',
            'aggv3_diagnosed_ndd',
            'aggv3_other',
            'deCODE_variant_carriers',
        ]
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for (gene_name, gene_id), stats in sorted(gene_stats.items(), key=lambda x: (x[0][0], x[0][1])):
            writer.writerow({
                'snoRNA_gene': gene_name,
                'gene_id': gene_id,
                'aggv3_undiagnosed_ndd': f"{len(stats['aggv3_undiagnosed'])}/{phenotype_denominators['undiagnosed NDD']}",
                'aggv3_diagnosed_ndd': f"{len(stats['aggv3_diagnosed'])}/{phenotype_denominators['diagnosed NDD']}",
                'aggv3_other': f"{len(stats['aggv3_other'])}/{phenotype_denominators['other']}",
                'deCODE_variant_carriers': len(stats['deCODE_participants']),
            })
    print('Done.', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Summarize snoRNA NDD status with GEL denominators and deCODE separation')
    parser.add_argument('--variants', required=True, help='Input variant TSV or TSV.gz file with ParticipantId, VariantId, chr, start, end, study')
    parser.add_argument('--phenotype', required=True, help='GEL phenotype TSV file with platekey, ndd, and case_solved')
    parser.add_argument('--gtf', required=True, help='GTF file path to identify snoRNA regions when variant annotation is missing')
    parser.add_argument('--out', required=True, help='Output TSV path')
    args = parser.parse_args()
    summarize_snoRNA_ndd_status(args.variants, args.phenotype, args.gtf, args.out)


if __name__ == '__main__':
    main()
