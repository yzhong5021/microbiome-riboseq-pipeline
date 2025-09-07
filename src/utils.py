"""
General helper functions for data processing-related tasks.

Functions:
    extract_species_name: given a file path, extracts and returns the species name.
    create_dir: creates a directory if it doesn't exist.
"""

import os

def extract_species_name(file_path):
    """given a file path, extracts and returns the species name."""

    basename = os.path.basename(file_path)

    return ' '.join(basename.split('_')[:-2]).strip()

def create_dir(out_path, dir_name):
    """creates a directory if it doesn't exist."""

    out_dir = os.path.join(out_path, dir_name)
    os.makedirs(out_dir, exist_ok=True)

    return out_dir
