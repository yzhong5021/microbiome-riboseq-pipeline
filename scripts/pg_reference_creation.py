#!/usr/bin/env python3

"""
creates reference centrifuger db made of pangenomes of all species > {thres} reads
"""

import argparse
import os
import subprocess
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

from src import cli_utils
from src import db_utils
from src import utils


def parse_args():
    """Parses command-line arguments."""
    return cli_utils.create_pg_reference_parser().parse_args()


def parse_taxonomy(taxo_path):
    """
    given taxonomy tsv with 1st col taxo id and 2nd col name,
    returns df mapping name to taxo id
    """
    taxo_df = pd.read_csv(taxo_path, sep='\t\t\t', header = None, names = ['taxo_id', 'name'])
    taxo_df = taxo_df.set_index('name')
    print("parsed taxonomy file. head:")
    print(taxo_df.head())
    return taxo_df


def create_file_mapping(pg_path, taxo_df, cluster_species_map, out_path):
    """
    creates the -l argument for centrifuger-build, a .tsv with first column reference path, second column taxo id
    """

    file_map_path = os.path.join(out_path, 'file_map.tsv')
    pg_path = os.path.join(out_path, 'pangenomes')

    with open(file_map_path, 'w') as f:
        for file in os.listdir(pg_path):
            clname = file.split('.')[0]
            filename = os.path.join(pg_path, file)
            taxo_id = taxo_df.loc[cluster_species_map[clname], 'taxo_id']
            f.write(f'{filename}\t{taxo_id}\n')

    return file_map_path


def run_cf(file_map_path, out_path):
    """
    creates the cf database.
    """
    nodes_dmp = r'/home/yochen/microbiome/cf/cf_taxonomy/nodes.dmp'
    names_dmp = r'/home/yochen/microbiome/cf/cf_taxonomy/names.dmp'
    cmd = [
        'centrifuger-build',
        '-l', file_map_path,
        '--taxonomy-tree', nodes_dmp,
        '--name-table', names_dmp,
        '--offrate', '3',
        '-o', os.path.join(out_path, 'pangenome'),
        '-t', '32']

    subprocess.run(cmd, check = True)


def main():
    args = parse_args()

    speci_df = pd.read_csv(args.speci_path, header=0, sep='\t')
    abun_df = pd.read_csv(args.agg_path, header=0, sep='\t')
    si_map = db_utils.get_specI_mapping(speci_df, abun_df, args.thres)
    
    cluster_species_map = {}
    for s in si_map:
        cluster_species_map[si_map[s][0]] = s

    species_list = list(si_map.keys())

    for species in tqdm(species_list, desc="Downloading species pangenomes..."):
        print(f"\n---{species}---")
        specI, _ = si_map[species]
        db_utils.download_pangenome_progenomes(specI, args.out_path)

    pg_dir = os.path.join(args.out_path, "pangenomes")

    taxo_df = parse_taxonomy(args.taxo_path)

    file_map_path = create_file_mapping(pg_dir, taxo_df, cluster_species_map, args.out_path)

    print("building centrifuge database...")

    run_cf(file_map_path, args.out_path)


if __name__ == "__main__":
    main()