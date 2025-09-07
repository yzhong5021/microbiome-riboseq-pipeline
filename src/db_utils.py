"""
Helper functions for creating, aligning, and manipulating reference sequence databases.

Functions:
    get_specI_mapping: Identifies species above a certain abundance threshold and maps them to the largest corresponding specI cluster from the progenomes database.
    download_pangenome: Downloads pangenome fasta and optionally eggNOG annotations for a specI cluster.
    download_pangenome_progenomes: Downloads pangenome from Progenomes database (legacy format).
    download_pangenome_with_size_check: Downloads pangenome and checks file size.
    build_bowtie2_index: Builds a bowtie2 index with prefix species name.
    align_reads: Aligns reads using bowtie2.
    process_sam_to_bam: Converts SAM to BAM, sorts, and indexes.
    get_bowtie2_stats: Gets bowtie2 alignment statistics from BAM file.
    get_mean_alignment_score: Calculates mean alignment score from BAM file.
"""

import os
import subprocess
from collections import defaultdict

import pandas as pd

from src import utils

def get_specI_mapping(speci_df, abun_df, thres):
    """
    Identifies species above a certain abundance threshold and maps them to the largest 
    corresponding specI cluster from the progenomes database.
    """
    cts = abun_df['count']
    nms = abun_df['name']

    species = set()
    for i in range(len(cts)):
        if cts[i] >= thres:
            species.add(nms[i])
    
    print(f"{len(species)} species passed threshold.")

    sp = speci_df['species']
    si = speci_df['specI']
    sz = speci_df['size.specI']
    res = defaultdict(lambda: ("", 0))
    for i in range(len(sp)):
        if sp[i] in species and sz[i] > res[sp[i]][1]:
            res[sp[i]] = (si[i], int(sz[i]))
            
    final_res = {}
    
    # correcting for taxonomic discrepancies (post-hoc)
    final_res['Stutzerimonas stutzeri'] = ('specI_v4_00181', 1)
    final_res['Sphingopyxis sp. FD7'] = ('specI_v4_36325', 1)
    final_res['Stutzerimonas kunmingensis'] = ('specI_v4_00449', 1)
    final_res['Ectopseudomonas oleovorans'] = ('specI_v4_01032', 1)
    final_res['Stutzerimonas chloritidismutans'] = ('specI_v4_00449', 1)
    final_res['Pseudomonas sp. OF001'] = ('specI_v4_08338', 1)
    final_res['Xenorhabdus szentirmaii DSM 16338'] = ('specI_v4_06794', 1)
    final_res['Pseudomonas marincola'] = ('specI_v4_02518', 1)
    final_res['Pseudomonas sp. 8Z'] = ('specI_v4_14874', 1)
    final_res['Pseudomonas aeruginosa'] = ('specI_v4_00059', 1)

    for name in species:
        if name in final_res:
            continue
        if res[name][1] > 0:
            final_res[name] = res[name]
        else:
            print(f'warning: no specI cluster for {name}')
    
    print(f"{len(final_res)} species passed to bowtie2.")
    return final_res

def download_pangenome(specI, pg_dir, include_eggnog=True):
    """
    Downloads pangenome fasta and optionally eggNOG annotations for a specI cluster.
    
    Args:
        specI: specI cluster identifier
        pg_dir: Directory to save files
        include_eggnog: Whether to download eggNOG annotations
        
    Returns:
        Tuple of (pangenome_fasta_path, eggnog_tsv_path) or (pangenome_fasta_path, None)
    """
    pangenome_fasta_path = os.path.join(pg_dir, f"{specI}.faa.gz")
    eggnog_tsv_path = os.path.join(pg_dir, f"{specI}.emapper.annotations.tsv") if include_eggnog else None
    
    pangenome_url = f"https://ftp.ebi.ac.uk/pub/databases/PGDB/latest/{specI}.faa.gz"
    
    if not os.path.exists(pangenome_fasta_path):
        print(f"Downloading pangenome for specI cluster {specI}...")
        subprocess.run(["wget", pangenome_url, "-O", pangenome_fasta_path], check=True, capture_output=True)

    if include_eggnog:
        eggnog_url = f"https://ftp.ebi.ac.uk/pub/databases/PGDB/latest/{specI}.emapper.annotations.tsv"
        if not os.path.exists(eggnog_tsv_path):
            print(f"Downloading eggNOG annotations for specI cluster {specI}...")
            subprocess.run(["wget", eggnog_url, "-O", eggnog_tsv_path], check=True, capture_output=True)

    return pangenome_fasta_path, eggnog_tsv_path


def download_pangenome_progenomes(specI, out_path):
    """
    Downloads pangenome from Progenomes database (legacy format).
    
    Args:
        specI: specI cluster identifier
        out_path: Output directory
        
    Returns:
        Path to downloaded pangenome directory
    """
    pg_dir = utils.create_dir(out_path, "progenomes")
    pangenome_url = f"https://progenomes.embl.de/dumpSequence.cgi?p={specI}&t=ga&a=specI"
    pangenome_fasta_path = os.path.join(pg_dir, f"{specI}.fa.gz")

    if not os.path.exists(pangenome_fasta_path):
        print(f"Downloading pangenome for specI cluster {specI}...")
        subprocess.run(["wget", pangenome_url, "-O", pangenome_fasta_path], check=True, capture_output=True)

    return pg_dir


def download_pangenome_with_size_check(specI, out_path):
    """
    Downloads pangenome and checks file size.
    
    Args:
        specI: specI cluster identifier
        out_path: Output directory
        
    Returns:
        Tuple of (pangenome_fasta_path, file_size_mb)
    """
    pg_dir = utils.create_dir(out_path, 'progenomes')
    pangenome_url = f"https://progenomes.embl.de/dumpSequence.cgi?p={specI}&t=ga&a=specI"
    pangenome_fasta_path = os.path.join(pg_dir, f"pangenome_{specI}.fa.gz")

    if not os.path.exists(pangenome_fasta_path):
        print(f"Downloading pangenome for specI cluster {specI}...")
        subprocess.run(["wget", pangenome_url, "-O", pangenome_fasta_path], check=True, capture_output=True)

    # Check file size
    cmd = f"gunzip -c {pangenome_fasta_path} | awk '!/^>/ {{s+=length($0)}} END {{print s/1000000}}'"
    mb_size_output = subprocess.run(
        cmd, 
        shell=True,
        check=True, 
        capture_output=True, 
        text=True
    )
    
    file_size_mb = float(mb_size_output.stdout.strip())
    print(f"Pangenome size: {file_size_mb:.2f} MB")
    
    return pangenome_fasta_path, file_size_mb

def build_bowtie2_index(prog_fasta, db_prefix):
    """
    Builds a bowtie2 index with prefix species name.
    """
    print(f"Building bowtie2 index...")
    subprocess.run(
        ["bowtie2-build", "-q", "-f", "--threads", "32", prog_fasta, db_prefix],
        check=True
    )

def align_reads(sub_fq, db_prefix, sam_file):
    """
    Aligns reads using bowtie2.
    """
    print(f"Aligning reads...")
    subprocess.run(
        ["bowtie2", "-U", sub_fq, "-x", db_prefix, "--very-sensitive", "-L", "15", "-N", "1", 
         "--no-unal", "-p", "32", "-S", sam_file],
        check=True
    )

def process_sam_to_bam(sam_file, bam_file, sorted_bam_file):
    """
    Converts SAM to BAM, sorts, and indexes.
    """
    print("Converting SAM to BAM, sorting, and indexing...")
    subprocess.run(["samtools", "view", "-b", sam_file, "-o", bam_file], check=True)
    os.remove(sam_file)
    subprocess.run(["samtools", "sort", bam_file, "-o", sorted_bam_file], check=True)
    os.remove(bam_file)
    subprocess.run(["samtools", "index", sorted_bam_file], check=True)

def get_bowtie2_stats(bam):
    """
    Gets bowtie2 alignment statistics from BAM file.
    """
    result = subprocess.run(
        ["samtools", "idxstats", "--threads", "8", bam],
        capture_output=True, text=True, check=True
    )
    return result.stdout

def get_mean_alignment_score(sorted_bam_path):
    """
    Calculates mean alignment score from BAM file.
    """
    cmd = (
        f"samtools view {sorted_bam_path} | "
        "awk '{for(i=12;i<=NF;i++) if($i~/^AS:i:/) {print substr($i,6); break}}' | "
        "awk '{s+=$1} END {if(NR>0) print s/NR; else print 0}'"
    )
    mean_as_proc = subprocess.run(
        cmd, shell=True, capture_output=True, text=True
    )
    return float(mean_as_proc.stdout.strip())

