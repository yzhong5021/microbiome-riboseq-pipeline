"""
Helper functions for FASTQ file processing.

Functions:
    get_subfq: Create subset FASTQ with all queries belonging to a specific species.
    extract_sequences_by_ids: Extract sequences from FASTQ file based on IDs in id_file.

"""

import os
import pandas as pd
from typing import Set

from src import utils


def get_subfq(query_map_df: pd.DataFrame, species: str, fastq_path: str, out_path: str) -> str:
    """Create subset FASTQ with all queries belonging to a specific species."""
    subfq_dir = utils.create_dir(out_path, "sub_fq")
    species_concat = "_".join(species.split())
    output_fastq = os.path.join(subfq_dir, f"{species_concat}.fq")
    
    queries = set(query_map_df[query_map_df['species_name'] == species]['query_id'])
    
    counter = 0
    with open(fastq_path, 'r') as f_in, open(output_fastq, 'w') as f_out:
        while True:
            header = f_in.readline()
            if not header:
                break
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()
            
            query_id = header.strip().split()[0][1:]
            if query_id in queries:
                counter += 1
                f_out.write(header)
                f_out.write(seq)
                f_out.write(plus)
                f_out.write(qual)
    
    print(f"{counter} sequences found in input fastq.")
    return output_fastq


def extract_sequences_by_ids(id_file: str, fastq_file: str, output_file: str) -> int:
    """Extract sequences from FASTQ file based on IDs in id_file."""

    with open(id_file, 'r') as f:
        target_ids = set(line.strip() for line in f)
    
    extracted_count = 0
    with open(fastq_file, 'r') as f_in, open(output_file, 'w') as f_out:
        while True:
            header = f_in.readline()
            if not header:
                break
            seq = f_in.readline()
            plus = f_in.readline()
            qual = f_in.readline()
            
            seq_id = header.strip().split()[0][1:]
            if seq_id in target_ids:
                f_out.write(header)
                f_out.write(seq)
                f_out.write(plus)
                f_out.write(qual)
                extracted_count += 1
    
    print(f"Extracted {extracted_count} sequences to {output_file}")
    return extracted_count
