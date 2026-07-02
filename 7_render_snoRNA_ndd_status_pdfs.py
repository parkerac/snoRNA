#!/usr/bin/env python3
"""
Render snoRNA NDD status tables as human-readable PDFs.

This script reads one or more TSV summary files and writes a PDF report for each.
It uses reportlab to create styled, wrapped table output suitable for review.

Usage:
  python 7_render_snoRNA_ndd_status_pdfs.py \
    --tsv gene_summary.tsv --out gene_summary.pdf \
    --tsv variant_summary.tsv --out variant_summary.pdf

If --out is omitted for a TSV file, the PDF name is inferred by replacing the
.tsv/.tsv.gz suffix with .pdf.
"""

import argparse
import csv
import os
import sys
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


PDF_PAGE_SIZE = landscape(A4)
MARGIN = 2 * cm


def infer_output_path(tsv_path):
    base = os.path.basename(tsv_path)
    if base.endswith('.tsv.gz'):
        base = base[:-7]
    elif base.endswith('.tsv'):
        base = base[:-4]
    return f'{base}.pdf'


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


def main():
    parser = argparse.ArgumentParser(description='Render snoRNA summary TSV files as human-readable PDFs')
    parser.add_argument('--tsv', required=True, action='append', help='Input TSV summary file')
    parser.add_argument('--out', action='append', help='Optional output PDF path; if omitted, inferred from TSV path')
    args = parser.parse_args()

    if args.out and len(args.out) != len(args.tsv):
        parser.error('If --out is provided, it must be given once per --tsv file')

    out_paths = args.out or [infer_output_path(tsv) for tsv in args.tsv]
    for tsv_path, out_path in zip(args.tsv, out_paths):
        build_pdf(tsv_path, out_path)


if __name__ == '__main__':
    main()
