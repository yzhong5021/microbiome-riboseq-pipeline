#!/bin/bash

kraken2-build --download-taxonomy --db k2_rs_db/ --threads 32

k2 add-to-library --db k2_rs_db --threads 32 --no-masking --file rf_bacteria/all_refseq_bacteria.fna

kraken2-build --download-library 'human' --db k2_rs_db/ --threads 64

kraken2-build --build --db k2_rs_db/ --kmer-len 25 --minimizer-len 23 --minimizer-spaces 2 --threads 64 --no-masking
