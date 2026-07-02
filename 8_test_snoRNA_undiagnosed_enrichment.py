#!/usr/bin/env python3
"""
Test for enrichment of snoRNA variant carriers in undiagnosed NDD cases.

This script reads a snoRNA variant TSV/TSV.gz file and a GEL phenotype TSV.
It counts carrier status by participant and performs Fisher exact tests for:
  - overall snoRNA carrier enrichment
  - per-gene enrichment
  - per-variant enrichment

The script uses exact inference to remain robust for sparse and small counts.
"""

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import defaultdict


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


def load_phenotypes(phenotype_path, platekey_col='platekey', ndd_col='ndd', case_solved_col='case_solved'):
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')
    phenotypes = {}
    totals = defaultdict(int)
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
            totals[group] += 1
    return phenotypes, totals


def fisher_exact(table):
    # table = [[a, b], [c, d]] with rows = [case, control] and cols = [carrier, noncarrier]
    a, b = table[0]
    c, d = table[1]
    if any(x < 0 for x in (a, b, c, d)):
        raise ValueError('Fisher exact input must be non-negative')
    total = a + b + c + d
    if total == 0:
        return float('nan'), 1.0, 1.0, 1.0

    def hypergeometric_p(x):
        return math.comb(a + b, x) * math.comb(c + d, a + c - x) / math.comb(total, a + c)

    observed = a
    row_sum = a + b
    col_sum = a + c
    low = max(0, row_sum - (c + d))
    high = min(row_sum, col_sum)
    observed_p = hypergeometric_p(observed)
    two_sided_p = 0.0
    less_p = 0.0
    greater_p = 0.0
    for x in range(low, high + 1):
        p = hypergeometric_p(x)
        if p <= observed_p + 1e-16:
            two_sided_p += p
        if x <= observed:
            less_p += p
        if x >= observed:
            greater_p += p

    if b == 0 or c == 0:
        odds_ratio = float('inf') if a * d > 0 else float('nan')
    else:
        odds_ratio = (a * d) / (b * c)
    return odds_ratio, two_sided_p, less_p, greater_p


def adjust_pvalues_bonferroni(p_values):
    n = len(p_values)
    return [min(p * n, 1.0) for p in p_values]


def adjust_pvalues_bh(p_values):
    n = len(p_values)
    sorted_indices = sorted(range(n), key=lambda i: p_values[i])
    adjusted = [0.0] * n
    prev_adj = 1.0
    for rank, idx in enumerate(sorted_indices, start=1):
        adj = p_values[idx] * n / rank
        prev_adj = min(prev_adj, adj)
        adjusted[idx] = min(prev_adj, 1.0)
    return adjusted


def adjust_pvalues(p_values, method):
    if method == 'none':
        return [min(p, 1.0) for p in p_values]
    if method == 'bonferroni':
        return adjust_pvalues_bonferroni(p_values)
    if method == 'bh':
        return adjust_pvalues_bh(p_values)
    raise ValueError(f'Unknown p-value adjustment method: {method}')


def build_carrier_sets(variants_path, phenotypes, exclude_studies=None):
    if exclude_studies is None:
        exclude_studies = set()
    participants_with_variants = defaultdict(lambda: {'genes': set(), 'variants': set()})

    with open_maybe_gzip(variants_path, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        required = {'ParticipantId', 'VariantId', 'study'}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Variants file must contain columns: {", ".join(sorted(missing))}')

        for row in reader:
            pid = row['ParticipantId']
            if pid not in phenotypes:
                continue
            study = str(row['study']).strip().lower()
            if study in exclude_studies:
                continue
            variant_id = row['VariantId']
            gene_name = row.get('snoRNA_gene') or row.get('snoRNA')
            gene_id = row.get('gene_id') or row.get('gene') or 'UNKNOWN'
            if gene_name:
                participants_with_variants[pid]['genes'].add((gene_name, gene_id))
                participants_with_variants[pid]['variants'].add((gene_name, gene_id, variant_id))
            else:
                participants_with_variants[pid]['variants'].add((None, None, variant_id))
    return participants_with_variants


def count_carriers(participant_ids, participants_with_variants, feature_type='gene'):
    counts = defaultdict(int)
    for pid in participant_ids:
        carrier_data = participants_with_variants.get(pid)
        if not carrier_data:
            continue
        if feature_type == 'gene':
            for gene_key in carrier_data['genes']:
                counts[gene_key] += 1
        elif feature_type == 'variant':
            for variant_key in carrier_data['variants']:
                counts[variant_key] += 1
        else:
            raise ValueError('feature_type must be gene or variant')
    return counts


def write_enrichment_results(out_path, rows, fieldnames):
    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def analyze_enrichment(variants_path, phenotype_path, out_path, case_group, control_group, exclude_studies, analysis_mode, min_carriers, p_adjust_method, alpha):
    phenotypes, phenotype_totals = load_phenotypes(phenotype_path)
    if case_group not in phenotype_totals:
        raise ValueError(f'Case group {case_group} not found in phenotype data')

    if control_group == 'other':
        control_groups = ['other']
    elif control_group == 'diagnosed NDD':
        control_groups = ['diagnosed NDD']
    elif control_group == 'all':
        control_groups = [group for group in phenotype_totals if group != case_group]
    else:
        control_groups = [control_group]

    case_participants = {pid for pid, group in phenotypes.items() if group == case_group}
    control_participants = {pid for pid, group in phenotypes.items() if group in control_groups}
    participants_with_variants = build_carrier_sets(variants_path, phenotypes, exclude_studies=exclude_studies)

    out_rows = []
    fieldnames = [
        'analysis_level',
        'snoRNA_gene',
        'gene_id',
        'VariantId',
        'case_group',
        'control_group',
        'case_participants',
        'control_participants',
        'case_carriers',
        'control_carriers',
        'case_carrier_freq',
        'control_carrier_freq',
        'case_noncarriers',
        'control_noncarriers',
        'odds_ratio',
        'p_value',
        'p_value_adjusted',
        'p_adjust_method',
        'fisher_p_two_sided',
        'fisher_p_less',
        'fisher_p_greater',
    ]

    total_case_carriers = sum(
        1 for pid in case_participants if participants_with_variants.get(pid, {}).get('variants')
    )
    total_control_carriers = sum(
        1 for pid in control_participants if participants_with_variants.get(pid, {}).get('variants')
    )
    case_noncarriers = len(case_participants) - total_case_carriers
    control_noncarriers = len(control_participants) - total_control_carriers
    overall_odds_ratio, overall_p_two_sided, overall_p_less, overall_p_greater = fisher_exact([
        [total_case_carriers, case_noncarriers],
        [total_control_carriers, control_noncarriers],
    ])
    out_rows.append({
        'analysis_level': 'overall',
        'snoRNA_gene': '',
        'gene_id': '',
        'VariantId': '',
        'case_group': case_group,
        'control_group': control_group,
        'case_participants': len(case_participants),
        'control_participants': len(control_participants),
        'case_carriers': total_case_carriers,
        'control_carriers': total_control_carriers,
        'case_carrier_freq': f'{total_case_carriers}/{len(case_participants)}',
        'control_carrier_freq': f'{total_control_carriers}/{len(control_participants)}',
        'case_noncarriers': case_noncarriers,
        'control_noncarriers': control_noncarriers,
        'odds_ratio': overall_odds_ratio,
        'p_value': overall_p_two_sided,
        'p_value_adjusted': overall_p_two_sided,
        'p_adjust_method': '',
        'fisher_p_two_sided': overall_p_two_sided,
        'fisher_p_less': overall_p_less,
        'fisher_p_greater': overall_p_greater,
    })

    if analysis_mode in {'gene', 'both'}:
        case_gene_counts = count_carriers(case_participants, participants_with_variants, feature_type='gene')
        control_gene_counts = count_carriers(control_participants, participants_with_variants, feature_type='gene')
        all_gene_keys = sorted(set(case_gene_counts) | set(control_gene_counts), key=lambda x: (x[0] or '', x[1] or ''))
        for gene_name, gene_id in all_gene_keys:
            case_carriers = case_gene_counts.get((gene_name, gene_id), 0)
            control_carriers = control_gene_counts.get((gene_name, gene_id), 0)
            if case_carriers + control_carriers < min_carriers:
                continue
            case_noncarriers = len(case_participants) - case_carriers
            control_noncarriers = len(control_participants) - control_carriers
            odds_ratio, p_two_sided, p_less, p_greater = fisher_exact([[case_carriers, case_noncarriers], [control_carriers, control_noncarriers]])
            p_raw = p_two_sided
            out_rows.append({
                'analysis_level': 'gene',
                'snoRNA_gene': gene_name,
                'gene_id': gene_id,
                'VariantId': '',
                'case_group': case_group,
                'control_group': control_group,
                'case_participants': len(case_participants),
                'control_participants': len(control_participants),
                'case_carriers': case_carriers,
                'control_carriers': control_carriers,
                'case_carrier_freq': f'{case_carriers}/{len(case_participants)}',
                'control_carrier_freq': f'{control_carriers}/{len(control_participants)}',
                'case_noncarriers': case_noncarriers,
                'control_noncarriers': control_noncarriers,
                'odds_ratio': odds_ratio,
                'p_value': p_raw,
                'p_value_adjusted': p_raw,
                'p_adjust_method': p_adjust_method,
                'fisher_p_two_sided': p_two_sided,
                'fisher_p_less': p_less,
                'fisher_p_greater': p_greater,
            })

    if analysis_mode in {'variant', 'both'}:
        case_variant_counts = count_carriers(case_participants, participants_with_variants, feature_type='variant')
        control_variant_counts = count_carriers(control_participants, participants_with_variants, feature_type='variant')
        all_variant_keys = sorted(
            set(case_variant_counts) | set(control_variant_counts),
            key=lambda x: (x[0] or '', x[1] or '', x[2] or ''),
        )
        for gene_name, gene_id, variant_id in all_variant_keys:
            case_carriers = case_variant_counts.get((gene_name, gene_id, variant_id), 0)
            control_carriers = control_variant_counts.get((gene_name, gene_id, variant_id), 0)
            if case_carriers + control_carriers < min_carriers:
                continue
            case_noncarriers = len(case_participants) - case_carriers
            control_noncarriers = len(control_participants) - control_carriers
            odds_ratio, p_two_sided, p_less, p_greater = fisher_exact([[case_carriers, case_noncarriers], [control_carriers, control_noncarriers]])
            p_raw = p_two_sided
            out_rows.append({
                'analysis_level': 'variant',
                'snoRNA_gene': gene_name,
                'gene_id': gene_id,
                'VariantId': variant_id,
                'case_group': case_group,
                'control_group': control_group,
                'case_participants': len(case_participants),
                'control_participants': len(control_participants),
                'case_carriers': case_carriers,
                'control_carriers': control_carriers,
                'case_carrier_freq': f'{case_carriers}/{len(case_participants)}',
                'control_carrier_freq': f'{control_carriers}/{len(control_participants)}',
                'case_noncarriers': case_noncarriers,
                'control_noncarriers': control_noncarriers,
                'odds_ratio': odds_ratio,
                'p_value': p_raw,
                'p_value_adjusted': p_raw,
                'p_adjust_method': p_adjust_method,
                'fisher_p_two_sided': p_two_sided,
                'fisher_p_less': p_less,
                'fisher_p_greater': p_greater,
            })

    if p_adjust_method not in {'none', 'bonferroni', 'bh'}:
        raise ValueError(f'Unknown p-value adjustment method: {p_adjust_method}')

    indices_to_adjust = [idx for idx, row in enumerate(out_rows) if row['analysis_level'] != 'overall']
    adjusted_pvalues = adjust_pvalues([out_rows[idx]['p_value'] for idx in indices_to_adjust], p_adjust_method)
    for idx, adj in zip(indices_to_adjust, adjusted_pvalues):
        out_rows[idx]['p_value_adjusted'] = adj

    filtered_rows = []
    for row in out_rows:
        threshold = row['p_value'] if row['analysis_level'] == 'overall' else row['p_value_adjusted']
        if threshold <= alpha:
            filtered_rows.append(row)

    write_enrichment_results(out_path, filtered_rows, fieldnames)


def main():
    parser = argparse.ArgumentParser(description='Test snoRNA variant enrichment in undiagnosed NDD cases')
    parser.add_argument('--variants', required=True, help='Input variant TSV or TSV.gz file with ParticipantId, VariantId, and snoRNA_gene or snoRNA')
    parser.add_argument('--phenotype', required=True, help='GEL phenotype TSV file with platekey, ndd, and case_solved')
    parser.add_argument('--out', required=True, help='Output TSV path for enrichment results')
    parser.add_argument('--case-group', default='undiagnosed NDD', help='Phenotype group to treat as cases')
    parser.add_argument('--control-group', default='other', choices=['other', 'diagnosed NDD', 'all'], help='Phenotype group to treat as controls')
    parser.add_argument('--exclude-studies', default='decode', help='Comma-separated study values to exclude from carrier counts')
    parser.add_argument('--analysis-mode', default='both', choices=['gene', 'variant', 'both'], help='Run gene-level, variant-level, or both enrichment tests')
    parser.add_argument('--min-carriers', type=int, default=1, help='Minimum total carriers for a feature to be reported')
    parser.add_argument('--p-adjust', default='bh', choices=['none', 'bonferroni', 'bh'], help='Multiple testing correction method')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance threshold after multiple testing correction')
    args = parser.parse_args()

    exclude_studies = {s.strip().lower() for s in args.exclude_studies.split(',') if s.strip()}
    analyze_enrichment(
        args.variants,
        args.phenotype,
        args.out,
        args.case_group,
        args.control_group,
        exclude_studies,
        args.analysis_mode,
        args.min_carriers,
        args.p_adjust,
        args.alpha,
    )


if __name__ == '__main__':
    main()
