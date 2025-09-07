#!/usr/bin/env python3
import sys

"""
This script collects a mapping of query IDs to species names from either a 
kraken2 or centrifuge alignment result file.

usage: python3 hit_cf_k2_ids.py <alignment_file> <species_taxids_file> <taxid_name_file> <output_file>
"""


if len(sys.argv) != 5:
    print("Usage: python3 script.py <input> <species_taxids> <taxid_to_name> <output_tsv>")
    sys.exit(1)

alignment_file, species_taxids_file, taxid_name_file, out_file = sys.argv[1:]

with open(species_taxids_file) as f:
    species_ids = { line.strip() for line in f if line.strip() }

taxid_to_name = {}
with open(taxid_name_file) as f:
    for line in f:
        parts = line.strip().split(None, 1)
        if len(parts) == 2:
            tid, name = parts
            taxid_to_name[tid] = name

found = 0
with open(alignment_file) as fin, open(out_file, 'w') as fout:
    fout.write("query_id\tspecies_name\n")
    next(fin)

    for line in fin:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 3:
            continue

        query_id = cols[0]
        taxid    = cols[2]

        if taxid in species_ids:
            name = taxid_to_name.get(taxid, "unknown")
            fout.write(f"{query_id}\t{name}\n")
            found += 1

print(f"wrote {found} records to {out_file}")