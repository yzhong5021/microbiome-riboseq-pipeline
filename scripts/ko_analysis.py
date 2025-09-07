#!/usr/bin/env python3

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root))

import argparse
import glob
import os
import shutil
import subprocess
from collections import defaultdict

import pandas as pd
from tqdm import tqdm

from src import cli_utils, db_utils, fastq_utils, func_utils, utils


def parse_args():
    """Parses command-line arguments."""
    return cli_utils.create_ko_analysis_parser().parse_args()


def bt2(prog_fasta, sub_fq, species, out_path):
    """
    builds a bt2 index with prefix species name and the following arguments:
        bowtie2-build -f --threads 32 {prog.fa} {out_path/bt2/db/{species name, concatenated by underscore}}
        stores the index in a directory out_path/bt2/db/
    then aligns sub_fq to the prog_fasta with the following arguments:
        bowtie2 -U {sub_fq} -x {bt2/idx_prefix} --very-sensitive -L 15 -N 1 --no-unal -p 32 -S out_path/bt2/res/{species name, conc. with underscore}.sam
    
    deletes the db directory afterward to save space
    uses samtools to compress output sam to bam
    returns the path to the .bam file
    """
    species_concat = "_".join(species.split())
    db_dir = utils.create_dir(out_path, os.path.join("bt2", "db", species_concat))
    res_dir = utils.create_dir(out_path, os.path.join("bt2", "res"))
    
    db_prefix = os.path.join(db_dir, species_concat)
    sam_file = os.path.join(res_dir, f"{species_concat}.sam")
    bam_file = os.path.join(res_dir, f"{species_concat}.bam")
    sorted_bam_file = os.path.join(res_dir, f"{species_concat}.sorted.bam")

    db_utils.build_bowtie2_index(prog_fasta, db_prefix)
    db_utils.align_reads(sub_fq, db_prefix, sam_file)
    
    shutil.rmtree(db_dir)
    db_utils.process_sam_to_bam(sam_file, bam_file, sorted_bam_file)
    
    return sorted_bam_file


def parse_bt2(bam, species, out_path):
    """
    uses samtools idxstats to get gene counts with the following command:
        samtools idxstats --threads 8 {bam}
    this will return a tsv (no header) with the following format:
        1st column = gene name; 2nd column = gene length; 3rd column = gene count; 4th column ignore
    normalize each count by gene length, normalizing to TPM; return a {species, concatenated}.tsv mapping gene to TPM to out_path/tpm
    also return another {species, concatenated}.tsv that maps gene to raw count/gene length to out_path/lnm_counts
    """
    species_concat = "_".join(species.split())
    tpm_dir = os.path.join(out_path, "tpm")
    lnm_dir = utils.create_dir(out_path, "lnm_counts")
    raw_dir = utils.create_dir(out_path, "raw_counts")
    tpm_dir = utils.create_dir(out_path, "tpm")

    result = db_utils.get_bowtie2_stats(bam)
    
    lines = result.strip().split('\n')
    data = [line.split('\t') for line in lines if not line.startswith('*')]
    df = pd.DataFrame(data, columns=['gene', 'length', 'count', 'unmapped'])
    df = df[['gene', 'length', 'count']].astype({'length': int, 'count': int})
    
    df = df[df['length'] > 0]
    
    df['rpk'] = df['count'] / (df['length'] / 1000)
    rpk_sum = df['rpk'].sum()
    df['tpm'] = (df['rpk'] / rpk_sum) * 1_000_000
    
    df['lnm_count'] = df['count'] / df['length']
    
    tpm_output = os.path.join(tpm_dir, f"{species_concat}.tsv")
    lnm_output = os.path.join(lnm_dir, f"{species_concat}.tsv")
    raw_output = os.path.join(raw_dir, f"{species_concat}.tsv")

    df[['gene', 'tpm']].to_csv(tpm_output, sep='\t', index=False)
    df[['gene', 'lnm_count']].to_csv(lnm_output, sep='\t', index=False)
    df[['gene', 'count']].to_csv(raw_output, sep='\t', index=False)
    
    return tpm_output, lnm_output, list(df['gene'])


def get_sets(kegg_tsv_path, species, out_path):
    """
    kegg_tsv_path specifies the path to a eggNOG annotation tsv path, with many columns but only two are important:
        QUERY_NAME and KEGG_KO.
    for every query, map to the KEGG_KO. write to out_path/kegg_sets/{species_concat}_kset.tsv

    repeat for pathway.
    """

    import re

    def filter_specific_pathways(pathway_string):
        if not isinstance(pathway_string, str):
            return 'Unknown'
        
        all_ids = re.split(r'[;,|]', pathway_string)
        
        valid_specific_pathways = []
        
        for pathway_id in all_ids:
            pathway_id = pathway_id.strip()
            match = re.search(r'(\d{5})', pathway_id)
            
            if match:
                path_num = int(match.group(1))
                if not (1100 <= path_num <= 1299):
                    valid_specific_pathways.append(pathway_id)
                    
        if valid_specific_pathways:
            return ','.join(valid_specific_pathways)
        else:
            return 'Unknown'

    species_concat = "_".join(species.split())

    ko_set_dir = utils.create_dir(out_path, "kegg_sets")
    ko_file = os.path.join(ko_set_dir, f"{species_concat}_kset.tsv")

    path_set_dir = utils.create_dir(out_path, "kpath_sets")
    path_file = os.path.join(path_set_dir, f"{species_concat}_pset.tsv")

    kegg_df = pd.read_csv(kegg_tsv_path, sep='\t', comment='#', header=0)

    ko_df = kegg_df[['QUERY_NAME', 'KEGG_KO']].copy()
    ko_df.fillna('Unknown', inplace = True)
    ko_df.rename(columns={'KEGG_KO': 'name'}, inplace = True)
    ko_df.to_csv(ko_file, sep='\t', index=False)
    
    path_df = kegg_df[['QUERY_NAME', 'KEGG_PATHWAY']].copy()
    path_df['KEGG_PATHWAY'] = path_df['KEGG_PATHWAY'].apply(filter_specific_pathways)
    path_df.rename(columns={'KEGG_PATHWAY': 'name'}, inplace = True)
    path_df.to_csv(path_file, sep='\t', index=False)

    return ko_file, path_file


def main():
    """
    1) for each species:
        a) get the speci clust identified by it
        b) download eggNOG KEGG terms/pangenome fasta from progenomes
        c) create subset fastq
        d) build/classify with bt2 index
        e) get length-normalized and TPM counts
        f) get KEGG set
        g) perform gsea
        h) delete temp files
    2) merge all counts/kegg terms
    3) gsea with merged counts/kegg terms
    with appropriate updates printed to stdout.
    """ 
    args = parse_args()
    
    utils.create_dir(args.out_path, "")

    speci_df = pd.read_csv(args.speci, header=0, sep='\t')
    abun_df = pd.read_csv(args.abun, header=0, sep='\t')
    query_map_df = pd.read_csv(args.q, header=0, sep='\t')
    
    si_map = db_utils.get_specI_mapping(speci_df, abun_df, args.thres)

    temp_df = pd.DataFrame.from_dict(si_map, orient='index', columns=['specI', 'size'])
    temp_df.index.name = 'species'
    temp_df.to_csv(os.path.join(args.out_path, "map_species_clust.tsv"), sep="\t")

    species_list = list(si_map.keys())

    for species in tqdm(species_list, desc="Processing species"):
        print(f"\nStarting processing for {species}...")
        species_concat = "_".join(species.split())

        ko_check = os.path.join(args.out_path, f"ko_gsea/{species_concat}_ko_gsea.tsv")
        path_check = os.path.join(args.out_path, f"pathway_gsea/{species_concat}_pathway_gsea.tsv")

        if os.path.exists(ko_check) and os.path.exists(path_check):
            print(f"ssgsea files found for '{species}'. Skipping.")
            continue

        # a
        specI, _ = si_map[species]

        # b
        pangenome_fasta, eggnog_tsv = db_utils.download_pangenome(specI, args.out_path)
        
        # c
        sub_fq = fastq_utils.get_subfq(query_map_df, species, args.f, args.out_path)

        # d
        bam_file = bt2(pangenome_fasta, sub_fq, species, args.out_path)

        # e
        tpm_tsv, _, gene_list = parse_bt2(bam_file, species, args.out_path)

        # f
        ko_set_file, path_set_file = get_sets(eggnog_tsv, species, args.out_path)

        # g
        if ko_set_file and os.path.exists(tpm_tsv):
            ko_path = os.path.join(args.out_path, 'ko_gsea')
            ko_gsea_path = func_utils.gsea(tpm_tsv, ko_set_file, 'ko', ko_path, species_concat)
            if ko_gsea_path:
                try:
                    func_utils.translate_ko(ko_gsea_path)
                except TypeError:
                    print(f"skipping {species}; no KO sets found.")
            else:
                print(f"skipping {species}; no KO sets found.")
        
        if path_set_file and os.path.exists(tpm_tsv):
            p_path = os.path.join(args.out_path, 'pathway_gsea')
            p_gsea_path = func_utils.gsea(tpm_tsv, path_set_file, 'pathway', p_path, species_concat)
            if p_gsea_path:
                try:
                    func_utils.translate_pathway(p_gsea_path)
                except TypeError:
                    print(f"skipping {species}; no pathway sets found.")
            else:
                print(f"skipping {species}; no pathway sets found.")

        print("Cleaning up intermediate files...")
        # h
        os.remove(pangenome_fasta)
        os.remove(eggnog_tsv)
        os.remove(sub_fq)
        os.remove(bam_file)
        os.remove(bam_file + ".bai")

    print("\n--- Merging all results ---")
    ko_merged_path = func_utils.merge_kocounts(
        os.path.join(args.out_path, "raw_counts"),
        os.path.join(args.out_path, "kegg_sets"),
        args.out_path
    )
    func_utils.translate_merged_ko(
        ko_merged_path
    )

    path_merged_path = func_utils.merge_pathwaycounts(
        os.path.join(args.out_path, "raw_counts"), 
        os.path.join(args.out_path, "kpath_sets"), 
        args.out_path
    )
    func_utils.translate_merged_pathway(
        path_merged_path
    )

    print("\nPipeline finished.")


if __name__ == "__main__":
    main()