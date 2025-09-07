#!/usr/bin/env python3

"""
Generate rarefaction curve data from a TSV file.

usage: python3 rarefaction_curve.py -f <filename> -r <rank> -i <index> -o <output_file>
"""

import argparse
import csv
import math
import multiprocessing as mp
import random
import sys
from collections import Counter


def shannon_diversity(counts):
    total = sum(counts)
    if total == 0:
        return 0.0

    diversity = 0.0
    for count in counts:
        if count > 0:
            p = count / total
            diversity -= p * math.log(p)
    return diversity


def simpson_index(counts):
    total = sum(counts)
    if total == 0 or total == 1:
        return 0.0

    sum_p_sq = 0.0
    for count in counts:
        sum_p_sq += (count * (count - 1))

    return 1.0 - (sum_p_sq / (total * (total - 1)))


def richness(counts):
    return len(counts)


def subsample(depth, population, weights):
    if depth == 0:
        return []
    choices = random.choices(population, weights=weights, k=depth)
    return list(Counter(choices).values())


def worker(args):
    d, population, weights, index_name = args
    counts = subsample(d, population, weights)

    if index_name == 'shannon':
        return shannon_diversity(counts)
    elif index_name == 'simpson':
        return simpson_index(counts)
    elif index_name == 'richness':
        return richness(counts)
    return 0.0


def rarefaction_curve(probs, total_counts, index_name, output_f):
    output_filename = output_f
    print(f"--> Output data will be written to: {output_filename}")

    means = []
    stds = []

    step = int(total_counts / 500) if total_counts > 20 else 1
    if step == 0: step = 1
    depths = range(0, total_counts + 1, step)

    population = list(range(len(probs)))
    num_reps = 500

    with mp.Pool(processes = 32) as pool:
        for i,d in enumerate(depths):
            if i % 100 == 0:
                print(str(int(i/100)) + "/5 complete.")
            worker_args = [(d, population, probs, index_name)] * num_reps
            reps = pool.map(worker, worker_args)

            n = len(reps)
            if n > 0:
                mean = sum(reps) / n
                sum_sq_diff = sum((x - mean) ** 2 for x in reps)
                std = math.sqrt(sum_sq_diff / n) if n > 1 else 0.0
                means.append(mean)
                stds.append(std)
            else:
                means.append(0.0)
                stds.append(0.0)

    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['depth', 'mean_diversity', 'std_dev'])
        for i, d in enumerate(depths):
            writer.writerow([d, means[i], stds[i]])

    return output_filename


def main():
    parser = argparse.ArgumentParser(
        description='Generate rarefaction curve data from a TSV file.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-f', '--filename', dest='filename', required=True, help='TSV file with a header and a count column')
    parser.add_argument('-r', '--rank', dest='rank', required=True, help='Rank (either species or genus)')
    parser.add_argument('-i', '--index', dest='index', help='Diversity index to calculate (shannon, simpson, richness)', default='shannon')
    parser.add_argument('-o', '--output', dest='output_f', help='Output location', required=True)
    args = parser.parse_args()

    index = args.index.lower()
    if index not in ['shannon', 'simpson', 'richness']:
        raise ValueError("Index must be one of: 'shannon', 'simpson', 'richness'")

    print(f"Reading counts from '{args.filename}'...")
    counts = []
    valid_count_headers = ['totalNumMatches', '  # Direct', 'count', 'new_est_reads']
    count_col = None
    
    with open(args.filename, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for col in valid_count_headers:
            if col in reader.fieldnames:
                count_col = col
                break

        for row in reader:
            counts.append(int(row[count_col]))

    total_counts = sum(counts)

    probs = [c / total_counts for c in counts]

    print(f"Total reads: {total_counts}")
    print(f"Number of {args.rank}: {len(counts)}")

    output_file = rarefaction_curve(probs, total_counts, index, args.output_f)

    print(f"\nRarefaction data generation complete; saved to {output_file}.")


if __name__ == "__main__":
    main()