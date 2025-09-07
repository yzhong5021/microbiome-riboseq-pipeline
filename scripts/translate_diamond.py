#!/usr/bin/env python3

"""
Aggregate and translate diamond search hits by designated taxonomic level.
Only the lowest-e-value hit per query (≤ threshold) is counted.
usage: python3 aggregate_taxonomy.py --map <map.sorted.tsv> --diamond <diamond.tsv> --level <species> --output <species_counts.tsv> --threads <8> --evalue <0.01>
"""

import argparse
import os
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate diamond search hits by taxonomic level. "
                    "Only the lowest-e-value hit per query (≤ threshold) is counted.",
    )
    parser.add_argument(
        "--map", required=True,
        help="Path to the mapping TSV file (protein_id<tab>organism_name)."
    )
    parser.add_argument(
        "--diamond", required=True,
        help="Path to the diamond search results TSV file."
    )
    parser.add_argument(
        "--output", required=True,
        help="Path for the output TSV file."
    )
    parser.add_argument(
        "--level", required=True, choices=["genus", "species"],
        help="The taxonomic level to aggregate counts on."
    )
    parser.add_argument(
        "--evalue", type=float, default=None,
        help="E-value threshold; only hits ≤ this are considered."
    )
    parser.add_argument(
        "--sort-by-count", action="store_true",
        help="If specified, sorts the output by count descending."
    )
    parser.add_argument(
        "--threads", type=int, default=os.cpu_count(),
        help="Number of threads for processing the map file."
    )
    return parser.parse_args()


def get_diamond_hit_counts(diamond_filepath, evalue_threshold):
    """
    For each query, keep only the hit with the lowest e-value (≤ threshold).
    """
    best_per_query = {}
    print(f"PASS 1: Scanning {diamond_filepath} for best hits...", file=sys.stderr)
    if evalue_threshold is not None:
        print(f"  using e-value cutoff <= {evalue_threshold}", file=sys.stderr)

    try:
        with open(diamond_filepath, 'r') as f:
            for i, line in enumerate(f, 1):
                parts = line.strip().split()
                if len(parts) < 3:
                    if line.strip():
                        print(f"Warning: malformed line {i}, skipping.", file=sys.stderr)
                    continue

                query_id, subject_acc, ev_str = parts[0], parts[1], parts[2]
                try:
                    ev = float(ev_str)
                except ValueError:
                    print(f"Warning: bad e-value on line {i}, skipping.", file=sys.stderr)
                    continue

                if evalue_threshold is not None and ev > evalue_threshold:
                    continue

                prev = best_per_query.get(query_id)
                if (prev is None) or (ev < prev[0]):
                    best_per_query[query_id] = (ev, subject_acc)
    except FileNotFoundError:
        print(f"Error: Diamond file not found at {diamond_filepath}", file=sys.stderr)
        sys.exit(1)

    subject_counts = defaultdict(int)
    for ev, subj in best_per_query.values():
        subject_counts[subj] += 1

    total_queries = len(best_per_query)
    unique_subjects = len(subject_counts)
    print(f"Selected best hit for {total_queries:,} queries across {unique_subjects:,} subjects.",
          file=sys.stderr)
    return subject_counts


def process_map_chunk(chunk, subject_counts, level):
    local_tax_counts = defaultdict(int)
    local_hits = 0

    for line in chunk:
        parts = line.rstrip("\n").split("\t", 1)
        if len(parts) != 2:
            continue
        subj_acc, org_name = parts
        cnt = subject_counts.get(subj_acc)
        if not cnt:
            continue

        toks = org_name.split()
        if level == 'genus':
            taxon = toks[0] if toks else org_name
        else:
            taxon = " ".join(toks[:2]) if len(toks) >= 2 else org_name

        local_tax_counts[taxon] += cnt
        local_hits += cnt

    return local_tax_counts, local_hits


def aggregate_from_map_multithreaded(map_filepath, subject_counts, level, num_threads, chunk_size=100000):
    print(f"\nPASS 2: Aggregating map file {map_filepath} with {num_threads} threads...", file=sys.stderr)
    final_counts = defaultdict(int)
    total_hits = 0

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        try:
            with open(map_filepath, 'r') as f:
                chunk = []
                for line in f:
                    chunk.append(line)
                    if len(chunk) >= chunk_size:
                        futures.append(executor.submit(process_map_chunk, chunk, subject_counts, level))
                        chunk = []
                if chunk:
                    futures.append(executor.submit(process_map_chunk, chunk, subject_counts, level))

            for fut in as_completed(futures):
                loc_counts, loc_hits = fut.result()
                total_hits += loc_hits
                for t, c in loc_counts.items():
                    final_counts[t] += c
        except FileNotFoundError:
            print(f"Error: Map file not found at {map_filepath}", file=sys.stderr)
            sys.exit(1)

    print(f"Aggregated {total_hits:,} total hits into {len(final_counts):,} unique taxa.", file=sys.stderr)
    return final_counts, total_hits


def write_output(output_filepath, counts, total_hits, sort_by_count):
    print(f"\nWriting output to {output_filepath}...", file=sys.stderr)
    if total_hits == 0:
        print("Warning: no hits to write—output will be empty.", file=sys.stderr)
        with open(output_filepath, 'w') as f:
            f.write("name\tcount\tproportion\n")
        return

    items = list(counts.items())
    if sort_by_count:
        items.sort(key=lambda x: x[1], reverse=True)

    with open(output_filepath, 'w') as f:
        f.write("name\tcount\tproportion\n")
        for name, cnt in items:
            prop = cnt / total_hits
            f.write(f"{name}\t{cnt}\t{prop:.6g}\n")
    print("Done.", file=sys.stderr)


def main():
    args = parse_args()
    subject_counts = get_diamond_hit_counts(args.diamond, args.evalue)
    if not subject_counts:
        print("No valid best hits found—exiting.", file=sys.stderr)
        write_output(args.output, {}, 0, args.sort_by_count)
        return

    final_counts, total_hits = aggregate_from_map_multithreaded(
        args.map, subject_counts, args.level, args.threads
    )
    write_output(args.output, final_counts, total_hits, args.sort_by_count)


if __name__ == "__main__":
    main()