"""
Command line interface utilities for metaribo-seq analysis.

Functions:
    create_base_parser: Create a base argument parser with common options.
    add_common_args: Add common arguments to a parser.
    add_species_args: Add species-specific arguments to a parser.
    add_sequence_args: Add sequence processing arguments to a parser.
    create_ko_analysis_parser: Create argument parser for KO analysis pipeline.
    create_specificity_analysis_parser: Create argument parser for specificity analysis pipeline.
    create_pg_reference_parser: Create argument parser for pangenome reference creation.
    create_aggregate_parser: Create argument parser for count aggregation.

"""

import argparse
from typing import Dict, Any, Optional


def create_base_parser(description: str) -> argparse.ArgumentParser:
    """Create a base argument parser with common options."""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    return parser


def add_common_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add common arguments to a parser."""
    parser.add_argument(
        "--output_dir", "-o", 
        default=".", 
        dest='out_path',
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--verbose", "-v", 
        action="store_true",
        help="Enable verbose output"
    )
    return parser


def add_species_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add species-specific arguments to a parser."""
    parser.add_argument(
        "--speci_map", "-s", 
        required=True, 
        dest='speci',
        help="Path to the progenomes specI mapping file"
    )
    parser.add_argument(
        "--abundance_tsv", "-a", 
        required=True, 
        dest='abun',
        help="Path to the abundance results file"
    )
    parser.add_argument(
        "--threshold", "-t", 
        default=5000, 
        type=int, 
        dest='thres',
        help="Threshold count for pangenome mapping (default: 5000)"
    )
    return parser


def add_sequence_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add sequence processing arguments to a parser."""
    parser.add_argument(
        "--fastq", "-f", 
        required=True, 
        dest='f',
        help="Path to raw FASTQ file"
    )
    parser.add_argument(
        "--query_mapping", "-q", 
        required=True, 
        dest='q',
        help="Path to species-to-query mapping file"
    )
    return parser


def create_ko_analysis_parser() -> argparse.ArgumentParser:
    """Create argument parser for KO analysis pipeline."""
    parser = create_base_parser(
        "Pangenome-based read mapping and gene set enrichment analysis pipeline"
    )
    parser = add_common_args(parser)
    parser = add_species_args(parser)
    parser = add_sequence_args(parser)
    return parser


def create_specificity_analysis_parser() -> argparse.ArgumentParser:
    """Create argument parser for specificity analysis pipeline."""
    parser = create_base_parser(
        "Pangenome-based read mapping and gene set enrichment analysis pipeline"
    )
    parser = add_common_args(parser)
    parser = add_species_args(parser)
    parser = add_sequence_args(parser)
    return parser


def create_pg_reference_parser() -> argparse.ArgumentParser:
    """Create argument parser for pangenome reference creation."""
    parser = create_base_parser(
        "Pipeline for building cf pangenome reference database"
    )
    parser = add_common_args(parser)
    parser.add_argument(
        "--taxonomy", "-s", 
        required=True, 
        dest='taxo_path',
        help="Path to the taxo id-species mapping file"
    )
    parser.add_argument(
        "--speci_map", "-p", 
        required=True, 
        dest='speci_path',
        help="Path to the progenomes specI mapping file"
    )
    parser.add_argument(
        "--abundance_tsv", "-a", 
        required=True, 
        dest='agg_path',
        help="Path to the abundance results file"
    )
    parser.add_argument(
        "--threshold", "-t", 
        default=5000, 
        type=int, 
        dest='thres',
        help="Threshold count for pangenome addition (default: 5000)"
    )
    return parser


def create_aggregate_parser() -> argparse.ArgumentParser:
    """Create argument parser for count aggregation."""
    parser = create_base_parser(
        "Aggregate gene counts to functional terms"
    )
    parser = add_common_args(parser)
    parser.add_argument(
        "--input_dir", "-i",
        required=True,
        help="Input directory containing raw counts and annotations"
    )
    parser.add_argument(
        "--ko_delimiter",
        default="[;,|]",
        help="Regex delimiter for KO terms (default: [;,|])"
    )
    parser.add_argument(
        "--pathway_delimiter", 
        default="[;,|]",
        help="Regex delimiter for pathway terms (default: [;,|])"
    )
    return parser
