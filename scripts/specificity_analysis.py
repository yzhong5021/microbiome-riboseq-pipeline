#!/usr/bin/env python3

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root))

import argparse
import os
import re
import shutil
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from src import cli_utils
from src import db_utils
from src import fastq_utils
from src import func_utils
from src import utils


def parse_args():
    """Parses command-line arguments."""
    return cli_utils.create_specificity_analysis_parser().parse_args()


def build_bt2(prog_fasta, species, out_path):
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

    db_dir = utils.create_dir(out_path, f'bt2/db/{species_concat}')
    res_dir = utils.create_dir(out_path, f'bt2/res/')

    db_prefix = os.path.join(db_dir, species_concat)

    print(f"Building bowtie2 index for {species}...")
    subprocess.run(
        ["bowtie2-build", "-q", "-f", "--threads", "32", prog_fasta, db_prefix],
        check=True
    )

    return db_prefix


def bt2(sub_fq, db_prefix, out_path, species = "", keep_bam = False):    
    log_path = os.path.join(out_path, r"bt2/temp.log")

    species_concat = "_".join(species.split())

    sam_path = os.path.join(out_path, r"bt2/" + species_concat + ".sam")

    with open(log_path, 'w') as f:
        subprocess.run(
            ["bowtie2", "-U", sub_fq, "-x", db_prefix, "--very-sensitive", "-L", "15", "-N", "1", 
            "--no-unal", "-p", "32", "-S", sam_path],
            stderr=f,
            check=True
        )
    
    try:
        log_content = open(log_path, 'r').read()
        print("bowtie2 log:", log_content)
        rate_match = re.search(r"(\d+\.\d+)% overall alignment rate", log_content)

        zero_match = re.search(r"^\s*(\d+).*aligned 0 times", log_content, re.MULTILINE)
        one_match = re.search(r"^\s*(\d+).*aligned exactly 1 time", log_content, re.MULTILINE)
        multi_match = re.search(r"^\s*(\d+).*aligned >1 times", log_content, re.MULTILINE)


        alignment_rate = float(rate_match.group(1)) / 100.0

        counts = np.array([
            int(zero_match.group(1)),
            int(one_match.group(1)),
            int(multi_match.group(1))
        ])

        alignment_score = db_utils.get_mean_alignment_score(sam_path)

        bam_path = sam_path.replace('.sam', '.bam')
        sorted_bam_path = sam_path.replace('.sam', '.sorted.bam')

        subprocess.run(["samtools", "view", "-b", sam_path, "-o", bam_path], check=True)
        subprocess.run(["samtools", "sort", bam_path, "-o", sorted_bam_path], check=True)
        subprocess.run(["samtools", "index", sorted_bam_path], check=True)

        if not keep_bam:
            os.remove(bam_path)
            os.remove(sorted_bam_path)
            os.remove(sorted_bam_path + ".bai")

        os.remove(sam_path)
        os.remove(log_path)

        return alignment_rate, counts, alignment_score

    except Exception as e:
        print(f"Error processing bowtie2 output: {e}")
        return 0.0, np.array([0, 0, 0]), 0.0


def main():
    """
    Main function for specificity analysis pipeline.
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

    alignment_results = defaultdict(dict)
    similarity_matrix = defaultdict(dict)
    score_matrix = defaultdict(dict)
    pg_sizes = {}

    subfq_paths = {}
    for species in species_list:
        sub_fq = fastq_utils.get_subfq(query_map_df, species, args.f, args.out_path)
        subfq_paths[species] = sub_fq

    species_agg_counts = defaultdict(lambda: [0] * 3)
    species_agg_ar = defaultdict(list)
    species_agg_arc = defaultdict(list)
    
    for species_pg in tqdm(species_list, desc="Aligning for species"):
        print(f"\n------- {species_pg} -------")
        specI, _ = si_map[species_pg]

        pangenome_fasta, pg_size = db_utils.download_pangenome_with_size_check(specI, args.out_path)
        pg_sizes[species_pg] = pg_size
        db_prefix = build_bt2(pangenome_fasta, species_pg, args.out_path)
        
        for species in species_list:
            # internal alignment
            if species == species_pg:
                alignment_rate, counts, alignment_score = bt2(
                    subfq_paths[species], db_prefix, args.out_path, species=species, keep_bam=True
                )
                alignment_results[species].update({
                    'int_al': alignment_rate,
                    'int_al_cor': alignment_rate / pg_size,
                    'int_as': alignment_score,
                    'int_unal_ct': counts[0],
                    'int_al1_ct': counts[1],
                    'int_al2_ct': counts[2],
                    'pg_size': pg_size
                })
            
            # external alignment
            else:
                alignment_rate, _, alignment_score = bt2(
                    subfq_paths[species], db_prefix, args.out_path, species=species, keep_bam=False
                )

            similarity_matrix[species][species_pg + '_pg'] = alignment_rate
            score_matrix[species][species_pg + '_pg'] = alignment_score

        os.remove(pangenome_fasta)
        shutil.rmtree(os.path.dirname(db_prefix))

    sm_df_temp = pd.DataFrame(similarity_matrix)
    for species in species_list:
        external_cols = [col for col in sm_df_temp.columns if not col.startswith(species)]
        external_rates = sm_df_temp.loc[species, external_cols].values
        
        external_rates_corrected = [
            rate / pg_sizes[col.replace('_pg', '')] for rate, col in zip(external_rates, external_cols)
        ]
        
        alignment_results[species].update({
            'ext_al_mu': np.mean(external_rates),
            'ext_al_std': np.std(external_rates),
            'ext_al_cor_mu': np.mean(external_rates_corrected),
            'ext_al_cor_std': np.std(external_rates_corrected)
        })

    ar_df = pd.DataFrame.from_dict(alignment_results, orient='index')
    sm_df = pd.DataFrame(similarity_matrix)
    as_df = pd.DataFrame(score_matrix)

    ar_path = os.path.join(args.out_path, "alignment_results.tsv")
    sm_path = os.path.join(args.out_path, "similarity_matrix.tsv")
    as_path = os.path.join(args.out_path, "score_matrix.tsv")

    ar_df.to_csv(ar_path, sep='\t')
    sm_df.to_csv(sm_path, sep='\t')
    as_df.to_csv(as_path, sep='\t')

    print(f"saved alignment results to {ar_path}, alignment similarity to {sm_path}, and score matrix to {as_path}.")


if __name__ == "__main__":
    main()