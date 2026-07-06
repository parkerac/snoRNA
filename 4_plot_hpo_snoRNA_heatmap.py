#!/usr/bin/env python3
"""
Plot HPO enrichment in small-RNA variant carriers versus non-carriers.

Input:
  - merged variant-phenotype TSV with ParticipantId and rna_gene/snoRNA columns
  - phenotype TSV for the full background cohort with platekey and hpo_terms

Output:
  - heatmap image (PNG) of log2 odds ratios
  - optional TSV matrix with log2 odds ratios
  - optional TSV with the underlying 2x2 counts

Usage:
  python 4_plot_hpo_snoRNA_heatmap.py \
    --merged snoRNA_variant_phenotype.tsv \
    --phenotype phenotype_data.tsv \
    --out-heatmap hpo_small_rna_enrichment_heatmap.png \
    --out-matrix hpo_small_rna_log2_or_matrix.tsv \
    --out-counts hpo_small_rna_enrichment_counts.tsv \
    --top 20

The heatmap value is log2(odds ratio) for each HPO term in carriers of a
small-RNA label versus all phenotyped participants who are not carriers of
that label. Positive values indicate enrichment among carriers.
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict

from snoRNA_utils import get_rna_class, get_rna_gene, require_columns


def split_hpo_terms(value, hpo_sep='auto'):
    if value is None:
        return set()
    value = str(value).strip()
    if not value:
        return set()
    if hpo_sep == 'auto':
        parts = re.split(r'[;,]', value)
    else:
        parts = value.split(hpo_sep)
    return {part.strip() for part in parts if part.strip()}


def make_rna_label(row):
    rna_gene = get_rna_gene(row)
    if not rna_gene:
        return None
    rna_class = get_rna_class(row)
    return f'{rna_class}:{rna_gene}' if rna_class else rna_gene


def load_carriers(merged_path):
    carriers_by_label = defaultdict(set)
    with open(merged_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {'ParticipantId'}, 'Merged file')
        if not {'rna_gene', 'snoRNA'} & set(reader.fieldnames or []):
            raise ValueError('Merged file must contain rna_gene or snoRNA column')

        for row in reader:
            label = make_rna_label(row)
            if label:
                carriers_by_label[label].add(row['ParticipantId'])
    return dict(carriers_by_label)


def load_phenotype_terms(phenotype_path, participant_col='platekey', hpo_col='hpo_terms', hpo_sep='auto'):
    participant_terms = {}
    with open(phenotype_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {participant_col, hpo_col}, 'Phenotype file')
        for row in reader:
            participant_id = row[participant_col]
            participant_terms.setdefault(participant_id, set()).update(split_hpo_terms(row.get(hpo_col), hpo_sep))
    return participant_terms


def log2_odds_ratio(a, b, c, d, pseudocount=0.5):
    odds_ratio = ((a + pseudocount) * (d + pseudocount)) / ((b + pseudocount) * (c + pseudocount))
    return math.log2(odds_ratio), odds_ratio


def build_enrichment(carriers_by_label, participant_terms, pseudocount=0.5):
    all_participants = set(participant_terms)
    all_terms = sorted({term for terms in participant_terms.values() for term in terms})
    matrix = {term: {} for term in all_terms}
    count_rows = []

    for label in sorted(carriers_by_label):
        carriers = carriers_by_label[label] & all_participants
        noncarriers = all_participants - carriers
        missing_carriers = carriers_by_label[label] - all_participants
        if missing_carriers:
            print(f'WARNING: {label} has {len(missing_carriers)} carriers absent from phenotype background', file=sys.stderr)

        for term in all_terms:
            a = sum(1 for participant_id in carriers if term in participant_terms[participant_id])
            b = len(carriers) - a
            c = sum(1 for participant_id in noncarriers if term in participant_terms[participant_id])
            d = len(noncarriers) - c

            if not carriers or not noncarriers:
                log2_or = float('nan')
                odds_ratio = float('nan')
            else:
                log2_or, odds_ratio = log2_odds_ratio(a, b, c, d, pseudocount=pseudocount)

            matrix[term][label] = log2_or
            count_rows.append({
                'rna_label': label,
                'hpo_term': term,
                'carriers_with_term': a,
                'carriers_without_term': b,
                'noncarriers_with_term': c,
                'noncarriers_without_term': d,
                'carrier_total': len(carriers),
                'noncarrier_total': len(noncarriers),
                'odds_ratio': odds_ratio,
                'log2_odds_ratio': log2_or,
            })

    return matrix, count_rows


def carrier_term_totals(count_rows):
    totals = defaultdict(int)
    for row in count_rows:
        totals[row['hpo_term']] += row['carriers_with_term']
    return totals


def select_top_terms(matrix, count_rows, top_n=20, top_by='carrier-count', min_carriers_with_term=1):
    carrier_totals = carrier_term_totals(count_rows)
    terms = [term for term in matrix if carrier_totals[term] >= min_carriers_with_term]

    if top_by == 'max-abs-log2-or':
        def score(term):
            values = [value for value in matrix[term].values() if not math.isnan(value)]
            return max((abs(value) for value in values), default=0.0)
    else:
        def score(term):
            return carrier_totals[term]

    return sorted(terms, key=lambda term: (-score(term), term))[:top_n]


def write_matrix(out_path, matrix, terms, labels):
    with open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['hpo_term', *labels])
        for term in terms:
            values = []
            for label in labels:
                value = matrix[term].get(label, float('nan'))
                values.append('' if math.isnan(value) else f'{value:.6g}')
            writer.writerow([term, *values])


def write_counts(out_path, count_rows, terms):
    term_set = set(terms)
    fieldnames = [
        'rna_label',
        'hpo_term',
        'carriers_with_term',
        'carriers_without_term',
        'noncarriers_with_term',
        'noncarriers_without_term',
        'carrier_total',
        'noncarrier_total',
        'odds_ratio',
        'log2_odds_ratio',
    ]
    with open(out_path, 'w', newline='') as out:
        writer = csv.DictWriter(out, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for row in count_rows:
            if row['hpo_term'] not in term_set:
                continue
            out_row = dict(row)
            for key in ('odds_ratio', 'log2_odds_ratio'):
                value = out_row[key]
                out_row[key] = '' if math.isnan(value) else f'{value:.6g}'
            writer.writerow(out_row)


def plot_heatmap(matrix, terms, labels, out_png, cmap='coolwarm'):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError('matplotlib is required to render the heatmap PNG') from exc

    data = [[matrix[term].get(label, float('nan')) for label in labels] for term in terms]
    finite_values = [value for row in data for value in row if not math.isnan(value)]
    if not finite_values:
        raise ValueError('No finite enrichment values to plot')

    bound = max(abs(min(finite_values)), abs(max(finite_values)), 1.0)
    fig_width = max(8, len(labels) * 0.45 + 3)
    fig_height = max(6, len(terms) * 0.35 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=200)
    image = ax.imshow(data, aspect='auto', cmap=cmap, vmin=-bound, vmax=bound)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms)
    ax.set_xlabel('small RNA')
    ax.set_ylabel('HPO term')
    fig.colorbar(image, ax=ax, label='log2 odds ratio')
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='Plot HPO enrichment in small-RNA variant carriers versus non-carriers')
    parser.add_argument('--merged', required=True, help='Merged TSV with rna_gene/snoRNA and ParticipantId columns')
    parser.add_argument('--phenotype', required=True, help='Full phenotype TSV with all background participants')
    parser.add_argument('--out-heatmap', required=True, help='Output PNG path for enrichment heatmap')
    parser.add_argument('--out-matrix', help='Optional TSV output of the log2 odds ratio matrix')
    parser.add_argument('--out-counts', help='Optional TSV output of 2x2 counts and odds ratios')
    parser.add_argument('--top', type=int, default=20, help='Number of HPO terms to display (default 20)')
    parser.add_argument('--top-by', choices=['carrier-count', 'max-abs-log2-or'], default='carrier-count', help='How to choose displayed HPO terms')
    parser.add_argument('--min-carriers-with-term', type=int, default=1, help='Minimum total carrier observations for an HPO term to be displayed')
    parser.add_argument('--participant-col', default='platekey', help='Participant ID column in phenotype file')
    parser.add_argument('--hpo-col', default='hpo_terms', help='HPO term column in phenotype file')
    parser.add_argument('--hpo-sep', default='auto', help='HPO separator; use "auto" for comma or semicolon')
    parser.add_argument('--pseudocount', type=float, default=0.5, help='Continuity correction for odds ratios')
    parser.add_argument('--cmap', default='coolwarm', help='Matplotlib colormap for the heatmap')
    args = parser.parse_args()

    for path, label in [(args.merged, 'Merged file'), (args.phenotype, 'Phenotype file')]:
        if not os.path.exists(path):
            print(f'{label} not found: {path}', file=sys.stderr)
            sys.exit(2)

    carriers_by_label = load_carriers(args.merged)
    if not carriers_by_label:
        print('No small-RNA carriers found in merged file.', file=sys.stderr)
        sys.exit(0)

    participant_terms = load_phenotype_terms(
        args.phenotype,
        participant_col=args.participant_col,
        hpo_col=args.hpo_col,
        hpo_sep=args.hpo_sep,
    )
    if not participant_terms:
        print('No participants found in phenotype file.', file=sys.stderr)
        sys.exit(0)

    matrix, count_rows = build_enrichment(carriers_by_label, participant_terms, pseudocount=args.pseudocount)
    labels = sorted(carriers_by_label)
    top_terms = select_top_terms(
        matrix,
        count_rows,
        top_n=args.top,
        top_by=args.top_by,
        min_carriers_with_term=args.min_carriers_with_term,
    )
    if not top_terms:
        print('No HPO terms passed the carrier-count filter.', file=sys.stderr)
        sys.exit(0)

    if args.out_matrix:
        write_matrix(args.out_matrix, matrix, top_terms, labels)
    if args.out_counts:
        write_counts(args.out_counts, count_rows, top_terms)

    plot_heatmap(matrix, top_terms, labels, args.out_heatmap, cmap=args.cmap)
    print(f'Heatmap saved to {args.out_heatmap}', file=sys.stderr)
    if args.out_matrix:
        print(f'Matrix saved to {args.out_matrix}', file=sys.stderr)
    if args.out_counts:
        print(f'Counts saved to {args.out_counts}', file=sys.stderr)


if __name__ == '__main__':
    main()
