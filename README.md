# aVirus_paper

This repository contains the analysis pipeline for viral detection and characterization in metagenomic data, with a focus on ancient DNA samples from coprolites. The repository is linked to the paper Denise et al 2025

## Repository Structure

- `02-scripts/`: Contains all analysis workflows and scripts
  - `PREP_fastp_merge_reads/`: Preprocessing of raw sequencing reads
  - `QUAL_aDNA_damage/`: Ancient DNA damage pattern analysis
  - `REFG_*_taxonomic_profile/`: Various taxonomic profiling tools (centrifuge, kraken2, krakenuniq, metabuli)
  - `REFG_binning_contigs/`: Metagenomic binning analysis
  - `REFG_genomad_detection/`: Viral genome detection using geNomad
  - `REFG_mapping_*/`: Read and contig mapping workflows
  - `PLOT_*/`: Visualization and figure generation scripts
  - `PLOT_notebook_figure_generation.ipynb`: Jupyter notebook for figure generation
  - `PLOT_beta_diversity.R`: R script for beta diversity analysis
- `03-data/`: Reference data and databases
- `04-analysis/`: Analysis results organized by data type
  - `empirical_data/`: Results from real sample analysis
  - `simulation/`: Results from simulated data analysis

## Key Features

- Multiple viral detection tools integration (geNomad, Kraken2, KrakenUniq, Centrifuge, Metabuli)
- Ancient DNA damage pattern analysis
- Read preprocessing and quality control
- Metagenomic assembly and binning
- Reference-based mapping
- Taxonomic profiling
- Visualization and result analysis

## Workflow Organization

Each analysis workflow in `02-scripts/` follows a standardized Snakemake structure:
- `config/`: Configuration files
- `workflow/`: 
  - `Snakefile`: Main workflow definition
  - `envs/`: Conda environment specifications
  - `rules/`: Individual workflow steps
  - `scripts/`: Custom analysis scripts
  - `schemas/`: JSON schemas for configuration validation


