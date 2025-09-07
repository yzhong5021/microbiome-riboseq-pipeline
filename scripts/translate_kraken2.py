#!/usr/bin/env python3
import argparse
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize Kraken2 report")

    p.add_argument("-i", "--report", required=True,
                   help="Kraken2 .report file")
    p.add_argument("-t", "--top", type=int, default=None,
                   help="Top N taxa to display (by # direct assignments)")
    p.add_argument("-o", "--output-file", required=True,
                   help="Path to write the summary (overwrites if exists)")
    p.add_argument("-r", "--rank", type=str, default = None,
                   help="Rank for summary. Use single-letter code (S for species, G for genus, etc.)")
    return p.parse_args()


def load_report(report_file):
    """
    Parse the Kraken2 report file. Returns a list of taxon entries and total reads.
    Each entry: {taxid, rank, name, covered, direct}
    """
    entries = []
    total_covered = None
    with open(report_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            try:
                pct_covered = float(parts[0])
                covered = int(parts[1])
                direct = int(parts[2])
                rank = parts[3]
                tid = int(parts[4])
                name = " ".join(parts[5:])
            except ValueError:
                continue
            # identify total at root
            if rank == 'R' and tid == 1:
                total_covered = covered
            entries.append({
                'taxid': tid,
                'rank': rank,
                'name': name,
                'covered': covered,
                'direct': direct
            })
    return entries, total_covered


def summarize(entries, total, top_n, rank):
    """
    Sort entries by direct count descending and take top_n. Compute proportions.
    Returns list of dicts with keys: name, rank, direct, covered, pct_direct, pct_covered
    """
    # sort
    sorted_entries = sorted(entries, key=lambda e: e['direct'], reverse=True)
    summary = []
    
    if top_n == None:
        for e in sorted_entries:
            if rank != None and e['rank'] != rank:
                continue
            direct = e['direct']
            covered = e['covered']
            pct_direct = (direct / total * 100) if total else 0
            pct_covered = (covered / total * 100) if total else 0
            summary.append({
                'name': e['name'],
                'rank': e['rank'],
                'direct': direct,
                'covered': covered,
                'pct_direct': pct_direct,
                'pct_covered': pct_covered
        })
    else:
        for e in sorted_entries[:top_n]:
            if rank != None and e['rank'] != rank:
                continue
            direct = e['direct']
            covered = e['covered']
            pct_direct = (direct / total * 100) if total else 0
            pct_covered = (covered / total * 100) if total else 0
            summary.append({
                'name': e['name'],
                'rank': e['rank'],
                'direct': direct,
                'covered': covered,
                'pct_direct': pct_direct,
                'pct_covered': pct_covered
            })
    return summary


def main():
    args = parse_args()
    entries, total = load_report(args.report)
    if total is None:
        raise RuntimeError("Could not determine total reads from report (missing root line)")
    summary = summarize(entries, total, args.top, args.rank)

    with open(args.output_file, 'w') as out:
        out.write(f"{'Scientific name':<30}\t{'Rank':<6}\t{'# Direct':>10}\t{'# Covered':>10}\t{'% Direct':>10}\t{'% Covered':>10}\n")
        for e in summary:
            out.write(f"{e['name']:<30}\t{e['rank']:<6}\t{e['direct']:>10}\t{e['covered']:>10}\t{e['pct_direct']:10.2f}\t{e['pct_covered']:10.2f}\n")


if __name__ == '__main__':
    main()