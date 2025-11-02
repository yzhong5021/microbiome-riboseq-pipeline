# Function Parameters and Commands

This file documents the key parameters used in the Metaribo-seq pilot study.  
These commands are not intended for replication, since sample data is not publicly available.

***

## Preprocessing
**Cutadapt**

-a AAAAAAAAAACAAAAAAAAAA --overlap=4 --trimmed-only
--maximum-length=56 --minimum-length=31 --quality-cutoff=28


## rRNA Removal
**SortMeRNA**

--ref SILVA_138.2_LSURef_NR99_tax_silva.fasta
--ref SILVA_138.2_SSURef_NR99_tax_silva.fasta
--reads input.fq.gz --aligned rrna_filtered.fq.gz
--other clean_input.fq.gz --fastx --sam --best
--num_alignments 1 --threads 32


## Classifiers
**Centrifuge**

-x reference -u clean_input.fq -k 1
--min-hitlen 15 -S classification_results.tsv -p 32


**Kraken2 build**

--kmer-len 25 --minimizer-len 22
--minimizer-spaces 3 --no-masking --threads 32

**DIAMOND blastx**

--db db.dmnd --query clean_input.fq
--out alignments.m8 --matrix BLOSUM62
--ultra-sensitive --masking 1 --top 1
--evalue 0.01 --threads 64

**Bowtie2**

-U <sub_fq> -x <db_prefix> --very-sensitive -L 15 -N 1 --no-unal -p 32 -S <sam_file>
