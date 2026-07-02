#!/usr/bin/env python3
"""Shared helpers for the snoRNA analysis scripts."""

import csv
import gzip
import os
import re
from collections import defaultdict


NDD_GROUPS = ('undiagnosed NDD', 'diagnosed NDD', 'other')


def open_maybe_gzip(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def require_columns(fieldnames, required, label='Input file'):
    missing = set(required) - set(fieldnames or [])
    if missing:
        raise ValueError(f'{label} must contain columns: {", ".join(sorted(missing))}')


def parse_gtf_attributes(attributes):
    parsed = {}
    for item in attributes.split(';'):
        item = item.strip()
        if not item or ' ' not in item:
            continue
        key, value = item.split(' ', 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def merge_intervals(intervals):
    intervals = sorted(intervals)
    if not intervals:
        return []

    merged = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
            continue
        merged.append((cur_start, cur_end))
        cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged


def parse_snoRNA_regions(gtf_path, feature_type='snoRNA'):
    """Return chrom -> sorted (start, end, gene_name, gene_id) snoRNA regions.

    Intervals are merged per gene identity, where identity is the pair
    (gene_name, gene_id). Grouping by gene_id avoids assigning a repeated
    snoRNA name's gene_id to variants at a different locus.
    """
    regions_by_gene = defaultdict(lambda: defaultdict(list))
    with open_maybe_gzip(gtf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue

            chrom, _, _, start, end, _, _, _, attributes = parts[:9]
            attr = parse_gtf_attributes(attributes)
            gene_type = attr.get('gene_type') or attr.get('gene_biotype')
            if gene_type != feature_type:
                continue

            gene_name = attr.get('gene_name', 'UNKNOWN')
            gene_id = attr.get('gene_id', 'UNKNOWN')
            regions_by_gene[(gene_name, gene_id)][chrom].append((int(start), int(end)))

    regions_by_chrom = defaultdict(list)
    for (gene_name, gene_id), chroms in regions_by_gene.items():
        for chrom, intervals in chroms.items():
            for start, end in merge_intervals(intervals):
                regions_by_chrom[chrom].append((start, end, gene_name, gene_id))

    return {chrom: sorted(regions) for chrom, regions in regions_by_chrom.items()}


def find_overlapping_snoRNAs(chrom, start, end, regions_by_chrom):
    """Return (gene_name, gene_id) pairs overlapping a 1-based inclusive interval."""
    if end < start:
        start, end = end, start

    matches = []
    for region_start, region_end, gene_name, gene_id in regions_by_chrom.get(chrom, []):
        if region_start > end:
            break
        if start <= region_end:
            matches.append((gene_name, gene_id))
    return matches


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
    if normalize_bool(ndd_value) is True:
        return 'diagnosed NDD' if normalize_case_solved(case_solved_value) == 'yes' else 'undiagnosed NDD'
    return 'other'


def load_ndd_phenotypes(phenotype_path, platekey_col='platekey', ndd_col='ndd', case_solved_col='case_solved'):
    if not os.path.exists(phenotype_path):
        raise FileNotFoundError(f'Phenotype file not found: {phenotype_path}')

    phenotypes = {}
    totals = defaultdict(set)
    with open(phenotype_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        require_columns(reader.fieldnames, {platekey_col, ndd_col, case_solved_col}, 'Phenotype file')
        for row in reader:
            participant_id = row[platekey_col]
            group = categorize_ndd_status(row.get(ndd_col), row.get(case_solved_col))
            phenotypes[participant_id] = group
            totals[group].add(participant_id)
    return phenotypes, {group: len(participants) for group, participants in totals.items()}


def ensure_ndd_denominators(denominators):
    for group in NDD_GROUPS:
        denominators.setdefault(group, 0)


def aggv3_status_key(group):
    if group == 'diagnosed NDD':
        return 'aggv3_diagnosed'
    if group == 'undiagnosed NDD':
        return 'aggv3_undiagnosed'
    return 'aggv3_other'


def natural_key(value):
    return tuple(int(part) if part.isdigit() else part for part in re.split(r'(\d+)', str(value)))
