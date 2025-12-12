# Analysis Scripts

This directory contains all the analysis workflows implemented as Snakemake pipelines. Each subdirectory represents a distinct analysis step in the viral metagenomic pipeline.

## Workflow Descriptions

### Data Preprocessing
- `PREP_fastp_merge_reads/`: Quality control and merging of paired-end reads
  - Uses fastp for adapter trimming and quality filtering
  - Merges paired-end reads when possible

### Quality Assessment
- `QUAL_aDNA_damage/`: Analysis of ancient DNA damage patterns
  - Uses tools like DamageProfiler
  - Includes ancient DNA authenticity assessment
  - Maps reads to reference sequences using BWA

### Reference-based Analysis (REFG)
- `REFG_centrifuge_taxonomic_profile/`: Taxonomic classification using Centrifuge
- `REFG_kraken2_taxonomic_profile/`: Taxonomic classification using Kraken2
- `REFG_krakenuniq_taxonomic_profile/`: Taxonomic profiling with KrakenUniq
- `REFG_metabuli_taxonomic_profile/`: Taxonomic classification using Metabuli
- `REFG_binning_contigs/`: Metagenomic binning workflow
  - Uses tools like MMseqs2, VAMB, and SemiBin
  - Includes contig size filtering and abundance estimation
- `REFG_concat_binning/`: Consolidation of binning results
  - Combines outputs from multiple binning approaches
  - Prepares unified bin sets for downstream annotation and quality control
- `REFG_genomad_detection/`: Viral genome detection
  - Uses geNomad for viral sequence identification
  - Includes CheckV for quality assessment of predicted viral sequences
- `REFG_genomad_simulation/`: Simulation and benchmarking for geNomad
  - Generates simulated viral/metagenomic read sets and reference perturbations
  - Evaluates geNomad sensitivity/specificity and creates benchmarking reports
- `REFG_mapping_contigs/`: Contig-based mapping analysis
  - Maps assembled contigs to reference databases for taxonomic and functional context
- `REFG_mapping_reads_bwa_sam2lca/`: Read-based mapping with taxonomic assignment
  - Uses BWA for read mapping
  - Assigns taxonomic labels to reads based on mapping results using sam2lca

### Visualization and Analysis
- `PLOT_beta_diversity.R`: R script to calculate beta diversity and generate visualizations
- `PLOT_contigs_over_reference/`: Visualization of contig coverage and alignment to reference genomes
- `PLOT_evaluate_simulated_data/`: Performance evaluation and plotting for simulated datasets

### Additional Resources
- `PLOT_notebook_figure_generation.ipynb`: Jupyter notebook for figure generation when interactive exploration is needed

Each workflow directory follows a standardized Snakemake layout and typically contains:
- `config/config.yaml`: Configuration settings and sample mappings
- `workflow/`:
  - `Snakefile`: Workflow definition
  - `envs/*.yaml`: Conda environment specifications for reproducibility
  - `rules/`: Modular workflow rules
  - `scripts/`: Custom Python/R scripts used by rules
  - `schemas/`: JSON/YAML schemas for configuration validation

When adding or running workflows, edit the corresponding `config/config.yaml` to point to your input data and reference databases.

