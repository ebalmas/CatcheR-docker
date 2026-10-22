#!/usr/bin/env python3
"""
plasmid_final_corrected.py
==========================
Drop-in replacement for catcheR's plasmid_final.sh + plasmid_final2.R
with the 4-nt offset bug fixed.

THE BUG (original plasmid_final.sh):
    cut -c 30-37  →  final_BC.txt   captures  GCGT + BARCODE[0:4]  (wrong)
    cut -c 38-43  →  final_UCI.txt  captures  BARCODE[4:8] + UCI[0:2]  (wrong)

THE FIX (this script):
    characters 34-41  →  full 8-nt BARCODE  (correct)
    characters 42-47  →  full 6-nt UCI      (correct)

Read structure (1-based character positions, forward strand):
    pos  1-13 : UMI   (13 nt)
    pos 14-29 : RED primer  CTGATCAGCGAGCTAC  (16 nt)
    pos 30-33 : GCGT  overhang  (4 nt)  ← old code started here (wrong)
    pos 34-41 : BARCODE  (8 nt)         ← correct start
    pos 42-47 : UCI      (6 nt)         ← correct
    pos 48+   : GGCGCGTTCATCTG...  anchor + tail

USAGE
-----
    python plasmid_final_corrected.py \\
        --fastq   your_file.fastq.gz \\
        --barcodes rc_barcodes_genes.csv \\
        --output  results/ \\
        [--DIs 1000]          minimum reads per clone for plots (default 1000)
        [--clones clones.txt] optional list of clones of interest (barcode_UCI format)
        [--umi_len   13]      UMI length in nt  (default 13)
        [--bc_start  34]      barcode start position, 1-based (default 34)
        [--bc_end    41]      barcode end position,   1-based (default 41)
        [--uci_start 42]      UCI start position,     1-based (default 42)
        [--uci_end   47]      UCI end position,       1-based (default 47)

INPUT FILES
-----------
    FASTQ:    plain .fastq or gzipped .fastq.gz
    Barcodes: CSV  barcode,name   (2-column, no header, comma-separated)
              e.g. ACTAGGAT,PTPN11_2

OUTPUTS (written to --output folder)
--------------------------------------
    final_BC.txt                   one barcode per read  (corrected, 8 nt)
    final_UCI.txt                  one UCI per read      (corrected, 6 nt)
    final_UMI.txt                  one UMI per read      (13 nt)
    complete_table.csv             all reads: UMI, barcode, UCI, clone, gene
    distribution_all_clones.csv    all clone counts  (one row per unique BC_UCI)
    percentages.csv                % of reads per clone above DIs threshold
    counts_per_barcode.tsv         reads + unique UCIs per barcode
    counts_per_gene.tsv            reads + unique UCIs per gene
    qc_summary.txt                 mapping statistics
    [if --clones provided]
    clones_of_interest.csv         counts for clones of interest

REQUIREMENTS
------------
    Python 3.6+  —  no external packages needed.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# FASTQ reader
# ---------------------------------------------------------------------------
def open_fastq(path):
    """Yield (header, sequence, quality) for every read."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        while True:
            header = fh.readline().rstrip()
            if not header:
                break
            seq  = fh.readline().rstrip()
            fh.readline()           # '+' line
            qual = fh.readline().rstrip()
            yield header, seq, qual


# ---------------------------------------------------------------------------
# Load reference barcodes
# ---------------------------------------------------------------------------
def load_barcodes(csv_path):
    """
    Read rc_barcodes_genes.csv.
    Accepts:
      2-column:  barcode, name             (original catcheR format)
      3-column:  barcode, name, revcomp    (your extended format — uses col 1)
    Header row detected and skipped automatically.
    Returns dict { barcode_sequence : shRNA_name }
    """
    bc_to_name = {}
    with open(csv_path, newline='') as f:
        for row in csv.reader(f):
            if not row:
                continue
            if row[0].strip().lower() in ('barcode', 'bc', 'sequence'):
                continue
            bc_to_name[row[0].strip()] = row[1].strip()
    return bc_to_name


# ---------------------------------------------------------------------------
# Main extraction
# ---------------------------------------------------------------------------
def extract(fastq_path, umi_len, bc_start, bc_end, uci_start, uci_end):
    """
    Extract UMI, BC, UCI from every read using fixed character positions.
    Positions are 1-based (like cut), converted to 0-based Python slices here.

    Returns list of dicts: {umi, bc, uci}
    """
    # convert 1-based inclusive to 0-based Python slice
    bc_sl  = slice(bc_start  - 1, bc_end)    # e.g. 34-41 → [33:41]
    uci_sl = slice(uci_start - 1, uci_end)   # e.g. 42-47 → [41:47]
    umi_sl = slice(0, umi_len)                # e.g. 1-13  → [0:13]

    records = []
    for _header, seq, _qual in open_fastq(fastq_path):
        records.append({
            'umi': seq[umi_sl],
            'bc' : seq[bc_sl],
            'uci': seq[uci_sl],
        })
    return records


# ---------------------------------------------------------------------------
# Write flat files  (mirrors catcheR's final_BC.txt / final_UCI.txt / final_UMI.txt)
# ---------------------------------------------------------------------------
def write_flat_files(records, output_dir):
    with open(os.path.join(output_dir, 'final_BC.txt'),  'w') as f_bc,  \
         open(os.path.join(output_dir, 'final_UCI.txt'), 'w') as f_uci, \
         open(os.path.join(output_dir, 'final_UMI.txt'), 'w') as f_umi:
        for r in records:
            f_bc.write(r['bc']  + '\n')
            f_uci.write(r['uci'] + '\n')
            f_umi.write(r['umi'] + '\n')


# ---------------------------------------------------------------------------
# Build complete table and compute counts
# ---------------------------------------------------------------------------
def build_tables(records, bc_to_name):
    """
    Returns:
        complete_table  list of dicts: umi, bc, uci, clone, gene
        clone_counts    dict  clone → count
        stats           dict
    """
    complete_table = []
    clone_counts   = defaultdict(int)
    stats          = defaultdict(int)

    for r in records:
        stats['total'] += 1
        bc  = r['bc']
        uci = r['uci']
        umi = r['umi']

        name = bc_to_name.get(bc)
        if name is None:
            stats['unmatched'] += 1
            gene = 'UNMATCHED'
            name = 'UNMATCHED'
        else:
            stats['matched'] += 1
            gene = name.rsplit('_', 1)[0]

        clone = f"{bc}_{uci}"
        clone_counts[clone] += 1
        complete_table.append({
            'UMI':    umi,
            'barcode': bc,
            'UCI':    uci,
            'clone':  clone,
            'name':   name,
            'gene':   gene,
        })

    return complete_table, clone_counts, stats


# ---------------------------------------------------------------------------
# Write outputs  (mirrors plasmid_final2.R)
# ---------------------------------------------------------------------------
def write_outputs(complete_table, clone_counts, bc_to_name, DIs,
                  clones_of_interest, output_dir):

    os.makedirs(output_dir, exist_ok=True)
    total_matched = sum(1 for r in complete_table if r['name'] != 'UNMATCHED')

    # ── complete_table.csv ───────────────────────────────────────────────────
    with open(os.path.join(output_dir, 'complete_table.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['UMI','barcode','UCI','clone','name','gene'])
        w.writeheader()
        w.writerows(complete_table)

    # ── distribution_all_clones.csv  (all clone counts, annotated) ──────────
    # Build clone → (barcode, uci, name, gene, count)
    clone_meta = {}
    for r in complete_table:
        c = r['clone']
        if c not in clone_meta:
            clone_meta[c] = {
                'clone':   c,
                'barcode': r['barcode'],
                'UCI':     r['UCI'],
                'name':    r['name'],
                'gene':    r['gene'],
                'Freq':    0,
            }
        clone_meta[c]['Freq'] += 1

    all_clones = sorted(clone_meta.values(), key=lambda x: -x['Freq'])

    with open(os.path.join(output_dir, 'distribution_all_clones.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['clone','barcode','UCI','name','gene','Freq'])
        w.writeheader()
        w.writerows(all_clones)

    # ── percentages.csv  (clones above DIs threshold, % per clone) ──────────
    above_DIs = [c for c in all_clones if c['Freq'] > DIs and c['name'] != 'UNMATCHED']
    total_above = sum(c['Freq'] for c in above_DIs)

    with open(os.path.join(output_dir, 'percentages.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['clone', 'barcode', 'UCI', 'name', 'gene', 'Freq', 'percentage'])
        for c in above_DIs:
            pct = c['Freq'] / total_above * 100 if total_above else 0
            w.writerow([c['clone'], c['barcode'], c['UCI'],
                        c['name'], c['gene'], c['Freq'], f"{pct:.4f}"])

    # ── counts_per_barcode.tsv ───────────────────────────────────────────────
    bc_reads = defaultdict(int)
    bc_ucis  = defaultdict(set)
    for c in all_clones:
        if c['name'] != 'UNMATCHED':
            bc_reads[c['name']] += c['Freq']
            bc_ucis[c['name']].add(c['UCI'])

    assigned = sum(bc_reads.values())
    all_names = sorted(set(bc_to_name.values()))

    with open(os.path.join(output_dir, 'counts_per_barcode.tsv'), 'w') as f:
        f.write("barcode\tgene\treads\tpct_of_assigned\tunique_UCIs\tnote\n")
        for name in all_names:
            gene = name.rsplit('_', 1)[0]
            n    = bc_reads.get(name, 0)
            u    = len(bc_ucis.get(name, set()))
            pct  = n / assigned * 100 if assigned else 0
            note = "NOT_DETECTED" if n == 0 else ""
            f.write(f"{name}\t{gene}\t{n}\t{pct:.3f}\t{u}\t{note}\n")

    # ── counts_per_gene.tsv ──────────────────────────────────────────────────
    gene_reads = defaultdict(int)
    gene_ucis  = defaultdict(set)
    for name, n in bc_reads.items():
        gene = name.rsplit('_', 1)[0]
        gene_reads[gene] += n
        gene_ucis[gene].update(bc_ucis[name])

    with open(os.path.join(output_dir, 'counts_per_gene.tsv'), 'w') as f:
        f.write("gene\treads\tpct_of_assigned\tunique_UCIs\n")
        for gene, n in sorted(gene_reads.items(), key=lambda x: -x[1]):
            pct = n / assigned * 100 if assigned else 0
            f.write(f"{gene}\t{n}\t{pct:.2f}\t{len(gene_ucis[gene])}\n")

    # ── clones_of_interest.csv  (optional) ──────────────────────────────────
    if clones_of_interest:
        ofint = [c for c in all_clones if c['clone'] in clones_of_interest]
        with open(os.path.join(output_dir, 'clones_of_interest.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['clone','barcode','UCI','name','gene','Freq'])
            w.writeheader()
            w.writerows(ofint)
        print(f"  Clones of interest found: {len(ofint)} / {len(clones_of_interest)}")

    return assigned, above_DIs, all_clones


# ---------------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------------
def write_qc(stats, bc_reads, bc_to_name, DIs, output_dir):
    total    = stats['total']
    assigned = stats['matched']
    n_ref    = len(set(bc_to_name.values()))
    n_seen   = sum(1 for n in bc_to_name.values() if bc_reads.get(n, 0) > 0)
    missing  = sorted(n for n in bc_to_name.values() if bc_reads.get(n, 0) == 0)

    lines = [
        "=== plasmid_final_corrected.py  QC summary ===",
        "",
        "Extraction positions (1-based, corrected from original plasmid_final.sh):",
        "  UMI     : characters  1–13  (cut -c  1-13)",
        "  BARCODE : characters 34–41  (cut -c 34-41)  ← was 30-37, now fixed",
        "  UCI     : characters 42–47  (cut -c 42-47)  ← was 38-43, now fixed",
        "",
        f"Total reads processed   : {total:>10,}",
        f"Matched to reference    : {assigned:>10,}  ({assigned/total*100:.1f}%)",
        f"Unmatched               : {stats['unmatched']:>10,}  ({stats['unmatched']/total*100:.1f}%)",
        "",
        f"Reference barcodes      : {n_ref:>10}",
        f"Barcodes detected       : {n_seen:>10}",
        f"Barcodes missing        : {n_ref - n_seen:>10}",
        f"DIs threshold used      : {DIs:>10,}",
    ]
    if missing:
        lines += ["", f"Missing barcodes ({len(missing)}):"] + [f"  {m}" for m in missing]
    else:
        lines += ["", "All reference barcodes detected ✓"]

    qc_path = os.path.join(output_dir, 'qc_summary.txt')
    with open(qc_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Corrected catcheR plasmid_final pipeline (fixes 4-nt BC offset).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--fastq',     required=True,
                        help='Input FASTQ or FASTQ.gz (read 1)')
    parser.add_argument('--barcodes',  required=True,
                        help='rc_barcodes_genes.csv  (2-column: barcode,name)')
    parser.add_argument('--output',    required=True,
                        help='Output folder (created if absent)')
    parser.add_argument('--DIs',       type=int, default=1000,
                        help='Min reads per clone for plots/percentages (default: 1000)')
    parser.add_argument('--clones',    default=None,
                        help='Optional .txt file: list of clones of interest (barcode_UCI, one per line)')

    # Positional parameters — change these if your construct ever changes
    parser.add_argument('--umi_len',   type=int, default=13,
                        help='UMI length in nt (default: 13)')
    parser.add_argument('--bc_start',  type=int, default=34,
                        help='Barcode start position, 1-based (default: 34)')
    parser.add_argument('--bc_end',    type=int, default=41,
                        help='Barcode end position,   1-based (default: 41)')
    parser.add_argument('--uci_start', type=int, default=42,
                        help='UCI start position,     1-based (default: 42)')
    parser.add_argument('--uci_end',   type=int, default=47,
                        help='UCI end position,       1-based (default: 47)')
    args = parser.parse_args()

    # ── Print settings ───────────────────────────────────────────────────────
    print("Settings:", file=sys.stderr)
    print(f"  UMI     : characters  1–{args.umi_len}", file=sys.stderr)
    print(f"  BARCODE : characters {args.bc_start}–{args.bc_end}", file=sys.stderr)
    print(f"  UCI     : characters {args.uci_start}–{args.uci_end}", file=sys.stderr)
    print(f"  DIs threshold : {args.DIs}", file=sys.stderr)
    print(file=sys.stderr)

    # ── Load reference barcodes ──────────────────────────────────────────────
    print(f"Loading barcodes from: {args.barcodes}", file=sys.stderr)
    bc_to_name = load_barcodes(args.barcodes)
    print(f"  {len(bc_to_name)} reference barcodes loaded.", file=sys.stderr)

    # ── Load clones of interest (optional) ───────────────────────────────────
    clones_of_interest = set()
    if args.clones:
        with open(args.clones) as f:
            clones_of_interest = {line.strip() for line in f if line.strip()}
        print(f"  {len(clones_of_interest)} clones of interest loaded.", file=sys.stderr)

    # ── Extract from FASTQ ───────────────────────────────────────────────────
    print(f"Processing: {args.fastq}", file=sys.stderr)
    records = extract(
        args.fastq,
        umi_len=args.umi_len,
        bc_start=args.bc_start,   bc_end=args.bc_end,
        uci_start=args.uci_start, uci_end=args.uci_end,
    )
    print(f"  {len(records):,} reads extracted.", file=sys.stderr)

    # ── Write flat files ─────────────────────────────────────────────────────
    os.makedirs(args.output, exist_ok=True)
    print(f"Writing flat files (final_BC.txt, final_UCI.txt, final_UMI.txt) ...", file=sys.stderr)
    write_flat_files(records, args.output)

    # ── Build tables ─────────────────────────────────────────────────────────
    print("Building count tables ...", file=sys.stderr)
    complete_table, clone_counts, stats = build_tables(records, bc_to_name)

    # ── Write analysis outputs ───────────────────────────────────────────────
    print(f"Writing analysis outputs to: {args.output}/", file=sys.stderr)
    assigned, above_DIs, all_clones = write_outputs(
        complete_table, clone_counts, bc_to_name,
        DIs=args.DIs,
        clones_of_interest=clones_of_interest,
        output_dir=args.output,
    )

    # ── QC summary ───────────────────────────────────────────────────────────
    bc_reads = defaultdict(int)
    for c in all_clones:
        if c['name'] != 'UNMATCHED':
            bc_reads[c['name']] += c['Freq']

    write_qc(stats, bc_reads, bc_to_name, args.DIs, args.output)
    print(f"\nDone. All results in: {args.output}/", file=sys.stderr)


if __name__ == '__main__':
    main()
