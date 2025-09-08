# Pilot Pipeline for Microbial Ribo-seq (Metaribo-seq)

This repository contains the code and notebooks for a pilot study exploring the feasibility of applying ribosome profiling (Ribo-seq) to microbial soil environments. This is the first application of Ribo-seq in this context and is currently based on a single pore-water microbiome sample.

For a presentation of major results, see [presentations/project_presentation.pdf](presentations/project_presentation.pdf).

***

## Objectives
1. Taxonomic classification of individual Ribo-seq reads.
2. Functional analysis of translationally active microbiome features.

***

## Workflow

### 1. Taxonomic Classification
- Preprocessing: adapter trimming with Cutadapt, rRNA removal with SortMeRNA (SILVA v138.2 SSU and LSU).
- Read classification with multiple classifiers:
  - Centrifuger (GTDB r226, NCBI nt, RefSeq HBAV)
  - Kraken2 (RefSeq HBAV)
  - DIAMOND (NCBI nr, UniRef100)
  - Parameters can be found in [METHODS.md](METHODS.md).
- Classifier comparison and aggregation.
- Pangenome-based refinement using the ProGenomes v3 database and Centrifuger.
- Species specificity validation via cross-pangenome alignment.

Scripts and notebooks:
- [notebooks/basic_visualization.ipynb](notebooks/basic_visualization.ipynb)
- [notebooks/aggregate_classification.ipynb](notebooks/aggregate_classification.ipynb)
- [scripts/pg_reference_creation.py](scripts/pg_reference_creation.py)
- [notebooks/final_classification.ipynb](notebooks/final_classification.ipynb)
- [scripts/rarefaction_curve.py](scripts/rarefaction_curve.py)
- [notebooks/rarefaction.ipynb](notebooks/rarefaction.ipynb)

---

### 2. Functional Analysis
- Align reads to species-specific pangenomes (Bowtie2).
- Annotate genes with EggNOG-Mapper.
- Map to KO terms and KEGG pathways.
- Functional scoring:
  - Aggregated KO/KEGG counts (identify population-level functional enrichment).
  - ssGSEA-based enrichment per species (GSVA in R) (identify species-level functional enrichment).
  - Community-level aggregated ssGSEA enrichment scores (identify common functional capabilities across species).

Relevant scripts and notebooks:
- [scripts/ko_analysis.py](scripts/ko_analysis.py)
- [scripts/agg_counts.py](scripts/agg_counts.py)
- [scripts/ko_ssgsea.r](scripts/ko_ssgsea.r)
- [notebooks/functional_analysis.ipynb](notebooks/functional_analysis.ipynb)

---

### 3. Other Analyses
- Saturation analysis using *UMI-tools*.
  - Visualization: [notebooks/sat_analysis.ipynb](notebooks/sat_analysis.ipynb).

---

## Limitations

- Ribo-seq is not suited for taxonomic abundance estimation.
- Single-sample proof of concept.
- Results are highly dependent on classifier and reference choices, including pangenome selection.
- Some species (e.g., *Xenorhabdus szentirmaii*) showed low specificity and were excluded from functional analysis.

## Next Steps

- Integration with metagenomics data
  - De novo assembly/annotation of MAGs for more sample-specific results
- Generation of replicates and temporal samples to perform differential analyses.
- Correlation with environmental and geochemical data.
- Benchmark against metatranscriptomics.
  - Multi-omic approach via translational efficiency calculations.

---

## Methods and Citations

Methods can be found in [METHODS.md](METHODS.md).

Citations can be found in [CITATIONS.md](CITATIONS.md).
