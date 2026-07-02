#!/usr/bin/env python3
"""
Render snoRNA NDD status tables as human-readable PDFs or PNGs.

This script reads one or more TSV summary files and writes a report image for each.
It uses reportlab for PDF output and matplotlib for PNG output.

Usage:
  python 7_render_snoRNA_ndd_status_pdfs.py \
    --tsv gene_summary.tsv --out gene_summary.png --format png \
    --tsv variant_summary.tsv --out variant_summary.png --format png

If --out is omitted for a TSV file, the output name is inferred by replacing the
.tsv/.tsv.gz suffix with .pdf or .png depending on --format.
"""

import argparse
import csv
import os
import sys
import textwrap

import matplotlib.pyplot as plt
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


PDF_PAGE_SIZE = landscape(A4)
MARGIN = 2 * cm


def infer_output_path(tsv_path, fmt='pdf'):
    base = os.path.basename(tsv_path)
    if base.endswith('.tsv.gz'):
        base = base[:-7]
    elif base.endswith('.tsv'):
        base = base[:-4]
    suffix = '.png' if fmt == 'png' else '.pdf'
    return f'{base}{suffix}'


def read_tsv(tsv_path):
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f'TSV file not found: {tsv_path}')
    with open(tsv_path, 'r', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError(f'No header found in {tsv_path}')
        rows = [row for row in reader]
    return reader.fieldnames, rows


def make_paragraph(text, style):
    if text is None:
        text = ''
    return Paragraph(str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'), style)


def wrap_text(text, width=20):
    if text is None:
        return ''
    return '\n'.join(textwrap.wrap(str(text), width=width))


def build_png(tsv_path, out_path):
    headers, rows = read_tsv(tsv_path)
    title = os.path.basename(tsv_path).replace('.tsv.gz', '').replace('.tsv', '')
    wrapped_headers = [wrap_text(col, width=18) for col in headers]
    cell_text = []
    for row in rows:
        cell_text.append([wrap_text(row.get(col, ''), width=18) for col in headers])

    num_rows = len(cell_text) + 1
    num_cols = len(headers)
    fig_width = max(10, num_cols * 1.5)
    fig_height = max(6, min(50, num_rows * 0.35 + 2))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=150)
    ax.axis('off')

    table = ax.table(
        cellText=cell_text,
        colLabels=wrapped_headers,
        cellLoc='left',
        loc='center',
        colColours=['#d9d9d9'] * num_cols,
        edges='closed',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)

    plt.title(f'Summary Report: {title}', fontsize=14, pad=20)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.02)
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote PNG: {out_path}', file=sys.stderr)


def build_pdf(tsv_path, out_path):
    headers, rows = read_tsv(tsv_path)
    title = os.path.basename(tsv_path).replace('.tsv.gz', '').replace('.tsv', '')
    doc = SimpleDocTemplate(
        out_path,
        pagesize=PDF_PAGE_SIZE,
        leftMargin=MARGIN,
        rightMargin=MARGIN,
        topMargin=MARGIN,
        bottomMargin=MARGIN,
    )

    title_style = ParagraphStyle(
        name='Title',
        fontSize=18,
        leading=22,
        spaceAfter=12,
        alignment=1,
    )
    header_style = ParagraphStyle(
        name='Header',
        fontSize=8,
        leading=10,
        alignment=1,
        spaceAfter=4,
    )
    cell_style = ParagraphStyle(
        name='Cell',
        fontSize=7,
        leading=9,
        alignment=0,
    )

    data = [[make_paragraph(col, header_style) for col in headers]]
    for row in rows:
        data.append([make_paragraph(row.get(col, ''), cell_style) for col in headers])

    page_width, page_height = PDF_PAGE_SIZE
    available_width = page_width - (MARGIN * 2)
    col_width = max(50, available_width / max(1, len(headers)))
    col_widths = [col_width] * len(headers)

    table = Table(data, colWidths=col_widths, repeatRows=1)
    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#d9d9d9')),
        ('GRID', (0, 0), (-1, -1), 0.25, colors.grey),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 3),
        ('RIGHTPADDING', (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
    ]))

    elements = [Paragraph(f'Summary Report: {title}', title_style), Spacer(1, 12), table]
    doc.build(elements)
    print(f'Wrote PDF: {out_path}', file=sys.stderr)


def build_report(tsv_path, out_path, fmt='pdf'):
    if fmt == 'png':
        build_png(tsv_path, out_path)
    else:
        build_pdf(tsv_path, out_path)


def main():
    parser = argparse.ArgumentParser(description='Render snoRNA summary TSV files as human-readable PDF or PNG reports')
    parser.add_argument('--tsv', required=True, action='append', help='Input TSV summary file')
    parser.add_argument('--out', action='append', help='Optional output path; if omitted, inferred from TSV path and format')
    parser.add_argument('--format', choices=['pdf', 'png'], default=None, help='Optional output format; will be inferred from output path if omitted')
    args = parser.parse_args()

    if args.out and len(args.out) != len(args.tsv):
        parser.error('If --out is provided, it must be given once per --tsv file')

    out_paths = []
    for tsv_path, out_path in zip(args.tsv, args.out or [None] * len(args.tsv)):
        if out_path:
            out_paths.append(out_path)
        else:
            fmt = args.format or 'pdf'
            out_paths.append(infer_output_path(tsv_path, fmt=fmt))

    for tsv_path, out_path in zip(args.tsv, out_paths):
        fmt = args.format or ('png' if out_path.lower().endswith('.png') else 'pdf')
        build_report(tsv_path, out_path, fmt=fmt)


if __name__ == '__main__':
    main()
