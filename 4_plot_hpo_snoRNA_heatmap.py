#!/usr/bin/env python3
"""
Plot a heatmap of HPO term vs small-RNA gene showing the number of unique participants.

Input:
  - merged variant-phenotype TSV with columns: rna_gene/snoRNA, VariantId, ParticipantId, hpo_terms
    where `hpo_terms` is a comma-separated list of HPO terms per participant.

Output:
  - heatmap image (PNG)
  - optional TSV matrix with counts (HPO_term x small-RNA label)

Usage:
  python 4_plot_hpo_snoRNA_heatmap.py \
    --merged snoRNA_variant_phenotype.tsv \
    --out-heatmap hpo_snoRNA_heatmap.png \
    --out-matrix hpo_snoRNA_matrix.tsv \
    --top 20

Dependencies:
  python3 -m pip install pandas seaborn matplotlib

The script counts unique participants per (HPO term, small-RNA label) pair, selects the top N HPO terms by total participants,
then plots a clustered or plain heatmap (unspecified clustering).
"""

import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_and_prepare(merged_path, hpo_sep=','):
    df = pd.read_csv(merged_path, sep='\t', dtype=str)
    required = {'ParticipantId', 'hpo_terms'}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"Merged file must contain columns: {required}. Found: {set(df.columns)}")
    if 'rna_gene' not in df.columns and 'snoRNA' not in df.columns:
        raise ValueError(f"Merged file must contain rna_gene or snoRNA column. Found: {set(df.columns)}")

    has_class = 'rna_class' in df.columns
    df = df.copy()
    if 'rna_gene' in df.columns:
        df['rna_gene'] = df['rna_gene'].fillna('')
        if 'snoRNA' in df.columns:
            df['rna_gene'] = df['rna_gene'].where(df['rna_gene'].str.strip() != '', df['snoRNA'].fillna(''))
    else:
        df['rna_gene'] = df['snoRNA'].fillna('')
    if has_class:
        df['rna_class'] = df['rna_class'].fillna('snoRNA')
        df['rna_label'] = df['rna_class'] + ':' + df['rna_gene']
    else:
        df['rna_label'] = df['rna_gene']

    df = df[['rna_label', 'ParticipantId', 'hpo_terms']].copy()
    df['hpo_terms'] = df['hpo_terms'].fillna('')

    # Split and explode
    df['hpo_terms'] = df['hpo_terms'].astype(str)
    df['hpo_list'] = df['hpo_terms'].str.split(hpo_sep)
    df = df.explode('hpo_list')
    # strip whitespace
    df['hpo_term'] = df['hpo_list'].astype(str).str.strip()
    df = df[df['hpo_term'] != '']

    # Drop duplicates of same participant-RNA-term (ensures unique participants)
    df = df.drop_duplicates(subset=['ParticipantId', 'rna_label', 'hpo_term'])
    return df


def build_matrix(df):
    # Count unique participants per (hpo_term, small-RNA label)
    counts = df.groupby(['hpo_term', 'rna_label'])['ParticipantId'].nunique().reset_index()
    matrix = counts.pivot(index='hpo_term', columns='rna_label', values='ParticipantId').fillna(0).astype(int)
    return matrix


def select_top_terms(matrix, top_n=20):
    term_sums = matrix.sum(axis=1)
    top_terms = term_sums.sort_values(ascending=False).head(top_n).index.tolist()
    return matrix.loc[top_terms]


def plot_heatmap(matrix, out_png, figsize=None, cmap='viridis'):
    if figsize is None:
        figsize = (max(6, matrix.shape[1] * 0.3 + 3), max(6, matrix.shape[0] * 0.3 + 3))
    plt.figure(figsize=figsize)
    sns.set(style='white')
    ax = sns.heatmap(matrix, cmap=cmap, linewidths=0.5, linecolor='gray')
    plt.xlabel('small RNA')
    plt.ylabel('HPO term')
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main():
    p = argparse.ArgumentParser(description='Plot HPO x small-RNA heatmap of unique participants')
    p.add_argument('--merged', required=True, help='Merged TSV with rna_gene/snoRNA, ParticipantId, hpo_terms')
    p.add_argument('--out-heatmap', required=True, help='Output PNG path for heatmap')
    p.add_argument('--out-matrix', required=False, help='Optional TSV output of the matrix')
    p.add_argument('--top', type=int, default=20, help='Number of top HPO terms to display (default 20)')
    p.add_argument('--hpo-sep', default=',', help='Separator for HPO terms in hpo_terms column (default: ,)')
    args = p.parse_args()

    if not os.path.exists(args.merged):
        print(f"Merged file not found: {args.merged}", file=sys.stderr)
        sys.exit(2)

    df = load_and_prepare(args.merged, hpo_sep=args.hpo_sep)
    if df.empty:
        print("No HPO terms found after parsing.", file=sys.stderr)
        sys.exit(0)

    matrix = build_matrix(df)
    if matrix.empty:
        print("Matrix is empty (no overlaps).", file=sys.stderr)
        sys.exit(0)

    top_matrix = select_top_terms(matrix, top_n=args.top)

    if args.out_matrix:
        top_matrix.to_csv(args.out_matrix, sep='\t')

    plot_heatmap(top_matrix, args.out_heatmap)
    print(f"Heatmap saved to {args.out_heatmap}", file=sys.stderr)
    if args.out_matrix:
        print(f"Matrix saved to {args.out_matrix}", file=sys.stderr)


if __name__ == '__main__':
    main()
