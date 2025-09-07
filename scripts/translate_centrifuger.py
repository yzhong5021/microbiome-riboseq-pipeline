#!/usr/bin/env python3

"""
translates and bins centrifuger counts into either genus or species count tsv files.

usage: python3 translate_centrifuger.py <classification.tsv> <nodes.dmp> <names.dmp> <rank>
"""

import csv
import sys
from collections import defaultdict


def parse_taxonomy_files(nodes_file, names_file):
    """
    Parses nodes.dmp and names.dmp to create mapping dictionaries.
    Returns three dictionaries: taxid -> rank, taxid -> name, and taxid -> parent.
    """
    taxid_to_rank = {}
    taxid_to_name = {}
    taxid_to_parent = {}

    print("--> Reading taxonomy files...")
    with open(nodes_file, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.split('|')]
            tax_id = parts[0]
            parent_id = parts[1]
            rank = parts[2]
            taxid_to_rank[tax_id] = rank
            taxid_to_parent[tax_id] = parent_id

    with open(names_file, 'r') as f:
        for line in f:
            if 'scientific name' in line:
                parts = [p.strip() for p in line.split('|')]
                tax_id = parts[0]
                name = parts[1]
                taxid_to_name[tax_id] = name
                
    return taxid_to_rank, taxid_to_name, taxid_to_parent


def aggregate_and_sort(input_file, taxid_to_rank, taxid_to_name, taxid_to_parent, rank):
    """
    Aggregates numMatches from a classification file. If rank is 'genus',
    it rolls up species counts into their parent genus. Then, it filters,
    sorts, and writes the results.
    """
    if rank not in ['species', 'genus']:
        print(f"Error: Rank must be 'species' or 'genus', not '{rank}'.", file=sys.stderr)
        sys.exit(1)

    counts = defaultdict(int)
    
    print(f"--> Reading and aggregating counts from '{input_file}'...")
    with open(input_file, 'r', newline='') as infile:
        reader = csv.reader(infile, delimiter='\t')
        try:
            # Skip header, or process first line if no header
            first_line = next(reader)
            if 'readID' not in first_line[0]:
                if len(first_line) >= 8 and first_line[2].strip() != '0':
                    counts[first_line[2].strip()] += int(first_line[7])
        except (StopIteration, IndexError):
            print("Warning: Input file might be empty or malformed.", file=sys.stderr)
        
        for row in reader:
            if len(row) >= 8 and row[2].strip() != '0':
                counts[row[2].strip()] += int(row[7])

    print(f"--> Processing aggregated data for rank '{rank}'...")
    output_data = []

    if rank == 'species':
        for tax_id, total_matches in counts.items():
            if taxid_to_rank.get(tax_id) == 'species':
                name = taxid_to_name.get(tax_id, f"Unknown_species_{tax_id}")
                output_data.append((name, total_matches))
    
    elif rank == 'genus':
        genus_counts = defaultdict(int)
        for tax_id, total_matches in counts.items():
            current_rank = taxid_to_rank.get(tax_id)
            if current_rank == 'genus':
                genus_counts[tax_id] += total_matches
            elif current_rank == 'species':
                parent_id = taxid_to_parent.get(tax_id)
                # Ensure the parent exists and is a genus before adding
                if parent_id and taxid_to_rank.get(parent_id) == 'genus':
                    genus_counts[parent_id] += total_matches
        
        for tax_id, total_matches in genus_counts.items():
            name = taxid_to_name.get(tax_id, f"Unknown_genus_{tax_id}")
            output_data.append((name, total_matches))

    # Sort the final data by count, descending
    sorted_data = sorted(output_data, key=lambda item: item[1], reverse=True)

    # Write the result to a new file
    output_file = f"sorted_{rank}_abundance.tsv"
    print(f"--> Writing sorted results to '{output_file}'...")
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['name', 'totalNumMatches'])
        writer.writerows(sorted_data)
        
    print(f"Successfully created sorted file: {output_file}")


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 script_name.py <classification.tsv> <nodes.dmp> <names.dmp> <rank>")
        print("Example: python3 script_name.py results.tsv cf_taxonomy/nodes.dmp cf_taxonomy/names.dmp species")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    nodes_filename = sys.argv[2]
    names_filename = sys.argv[3]
    target_rank = sys.argv[4].lower()
    
    try:
        rank_map, name_map, parent_map = parse_taxonomy_files(nodes_filename, names_filename)
        aggregate_and_sort(input_filename, rank_map, name_map, parent_map, target_rank)
    except FileNotFoundError as e:
        print(f"Error: A required file was not found.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)