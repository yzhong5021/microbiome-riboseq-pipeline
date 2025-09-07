#!/usr/bin/env python3

"""
Collects a mapping of query IDs to species names from a DIAMOND alignment TSV,
given an accession to species name mapping file and an e-value cutoff.

usage: python3 hit_diamond_ids.py <diamond.tsv> <acc_to_species.tsv> <output.tsv> <eval_cutoff>
"""

import sys

def pick_best_hits(diamond_path, eval_cutoff):
    """
    Pass 1: for each queryID keep only the hit with the smallest e-value
    less than or equal to the evalue cutoff.
    Returns dict: queryID -> best_accID, and set of all best_accIDs.
    """
    best = {}      # queryID -> (best_ev, accID)
    with open(diamond_path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            qid, acc, ev_s = cols[0], cols[1], cols[2]
            try:
                ev = float(ev_s)
            except ValueError:
                continue
            if ev > eval_cutoff:
                continue

            prev = best.get(qid)
            if prev is None or ev < prev[0]:
                best[qid] = (ev, acc)

    query2acc = {}
    needed_accs = set()
    for qid, (ev, acc) in best.items():
        query2acc[qid] = acc
        needed_accs.add(acc)

    return query2acc, needed_accs


def load_species_names(acc2spec_path, needed_accs):
    """
    Stream acc_to_species.tsv for needed accessions to save memory.
    """
    acc2name = {}
    needs = set(needed_accs)
    with open(acc2spec_path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t", 1)
            if len(cols) != 2:
                continue
            acc, name = cols
            if acc in needs:
                if len(name.split()) >= 2:
                    acc2name[acc] = name
                needs.remove(acc)
                if not needs:
                    break
    return acc2name


def write_output(output_path, query2acc, acc2name):
    with open(output_path, "w") as out:
        for qid, acc in query2acc.items():
            name = acc2name.get(acc)
            if name:
                out.write(f"{qid}\t{name}\n")


def main():
    if len(sys.argv) != 5:
        print(f"Usage: python3 hit_diamond_ids.py "
              "<diamond.tsv> <acc_to_species.tsv> <output.tsv> <eval_cutoff>", file=sys.stderr)
        sys.exit(1)

    diamond_tsv    = sys.argv[1]
    acc2spec_tsv   = sys.argv[2]
    output_tsv     = sys.argv[3]
    eval_cutoff    = float(sys.argv[4])

    query2acc, needed_accs = pick_best_hits(diamond_tsv, eval_cutoff)
    print(f"Found best hits for {len(query2acc):,} queries, "
          f"{len(needed_accs):,} unique accIDs.", file=sys.stderr)

    if not needed_accs:
        print("no hits passed the e-value cutoff.", file=sys.stderr)
        open(output_tsv, "w").close()
        return

    acc2name = load_species_names(acc2spec_tsv, needed_accs)
    print(f"Mapped {len(acc2name):,} of {len(needed_accs):,} accIDs to species names.", file=sys.stderr)

    write_output(output_tsv, query2acc, acc2name)
    print(f"Wrote {len(acc2name):,} lines to {output_tsv}.", file=sys.stderr)


if __name__ == "__main__":
    main()