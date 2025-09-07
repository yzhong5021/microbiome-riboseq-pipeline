#!/usr/bin/env python3

"""
Aggregates raw gene counts to functional terms (either KEGG KOs or Pathways) for individual species.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root))

import glob
import os

import pandas as pd

from src import cli_utils
from src import func_utils
from src import utils


def aggregate_ko_counts(raw_counts_file, ko_annotation_file, out_path, species_concat, ko_delim='[;,|]'):
    """
    Aggregates raw gene counts to KEGG Orthology (KO) counts for a single species.
    """

    out_dir = utils.create_dir(out_path, "raw_spcts_ko")

    counts_df = pd.read_csv(raw_counts_file, sep='\t', header=0)
    ko_df = pd.read_csv(ko_annotation_file, sep='\t', header=0, usecols=['QUERY_NAME', 'name'])

    merged = counts_df.merge(ko_df, left_on="gene", right_on="QUERY_NAME", how="inner")

    merged.rename(columns={'name': 'ko_id'}, inplace=True)
    merged["ko_id"] = merged["ko_id"].astype(str).str.split(ko_delim, regex=True)
    exploded = merged.explode("ko_id")
    exploded["ko_id"] = exploded["ko_id"].str.strip()

    exploded = exploded[exploded['ko_id'].str.lower().notna() & (exploded['ko_id'].str.lower() != 'unknown') & (exploded['ko_id'] != '')]

    ko_counts = (
        exploded.groupby("ko_id", as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )

    output_file = os.path.join(out_dir, f"{species_concat}_ko_rc.tsv")
    ko_counts.to_csv(output_file, sep='\t', index=False)

    print(f"Aggregated {len(ko_counts)} KOs for {species_concat} -> {output_file}")
    return output_file


def aggregate_pathway_counts(raw_counts_file, pathway_annotation_file, out_path, species_concat, path_delim='[;,|]'):
    """
    Aggregates raw gene counts to KEGG Pathway (KP) counts for a single species.
    """
    
    out_dir = utils.create_dir(out_path, 'raw_spcts_kp')

    counts_df = pd.read_csv(raw_counts_file, sep='\t', header=0)
    path_df = pd.read_csv(pathway_annotation_file, sep='\t', header=0, usecols=['QUERY_NAME', 'name'])

    merged = counts_df.merge(path_df, left_on="gene", right_on="QUERY_NAME", how="inner")

    merged.rename(columns={'name': 'pathway_id'}, inplace=True)
    merged["pathway_id"] = merged["pathway_id"].astype(str).str.split(path_delim, regex=True)
    exploded = merged.explode("pathway_id")

    exploded = exploded[exploded['pathway_id'].str.lower().notna() & (exploded['pathway_id'].str.lower() != 'unknown') & (exploded['pathway_id'] != '')]
    exploded["pathway_id"] = exploded["pathway_id"].str.strip()

    exploded = exploded[~exploded['pathway_id'].str.contains('ko', case=False, na=False)]

    pathway_counts = (
        exploded.groupby("pathway_id", as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )

    output_file = os.path.join(out_dir, f"{species_concat}_kp_rc.tsv")
    pathway_counts.to_csv(output_file, sep='\t', index=False)

    print(f"Aggregated {len(pathway_counts)} pathways for {species_concat} -> {output_file}")
    return output_file


def main():

    base_path = os.path.dirname(os.path.abspath(__file__))
    counts_dir = os.path.join(base_path, "raw_counts")
    kegg_sets_dir = os.path.join(base_path, "kegg_sets")
    kpath_sets_dir = os.path.join(base_path, "kpath_sets")
    output_path = base_path

    if not os.path.isdir(counts_dir):
        print(f"Error: Input directory '{counts_dir}' not found.")
        return

    count_files = glob.glob(os.path.join(counts_dir, "*.tsv"))

    print(f"Found {len(count_files)} species files to process.")

    for raw_counts_file in count_files:
        species_name = utils.extract_species_name(raw_counts_file)
        print(f"\n--- Processing: {species_name} ---")

        ko_annotation_file = os.path.join(kegg_sets_dir, f"{species_name}_kset.tsv")
        pathway_annotation_file = os.path.join(kpath_sets_dir, f"{species_name}_pset.tsv")

        if os.path.exists(ko_annotation_file):
            ko_output_file = aggregate_ko_counts(raw_counts_file, ko_annotation_file, output_path, species_name)
            func_utils.translate_ko_ids(ko_output_file)
        else:
            print(f"Warning: KO annotation file not found for {species_name}. Skipping KO aggregation.")
            print(f"  (Looked for: {ko_annotation_file})")

        if os.path.exists(pathway_annotation_file):
            pathway_output_file = aggregate_pathway_counts(raw_counts_file, pathway_annotation_file, output_path, species_name)
            func_utils.translate_pathway_ids(pathway_output_file)
        else:
            print(f"Warning: Pathway annotation file not found for {species_name}. Skipping pathway aggregation.")
            print(f"  (Looked for: {pathway_annotation_file})")

    print("\nScript finished.")


if __name__ == '__main__':
    main()