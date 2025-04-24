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
  - Uses tools like MMseqs2 and VAMB and SemiBin
  - Includes contig size filtering and abundance estimation
- `REFG_genomad_detection/`: Viral genome detection
  - Uses geNomad for viral sequence identification
  - Includes CheckV for quality assessment
- `REFG_mapping_contigs/`: Contig-based mapping analysis using vclust
- `REFG_mapping_reads_bwa_sam2lca/`: Read-based mapping with taxonomic assignment 
  - Uses BWA for read mapping
  - Assigns taxonomic labels to reads based on mapping results using sam2lca

### Visualization and Analysis
- `PLOT_beta_diversity.R`: Rscript to calculate Beta diversity analysis and visualization
- `PLOT_contigs_over_reference/`: Visualization of contig coverage
- `PLOT_evaluate_simulated_data/`: Performance evaluation using simulated data

### Additional Resources
- `PLOT_notebook_figure_generation.ipynb`: Jupyter notebook for figure generation if not done by external snakemake pipeline or Rscript

Each workflow directory contains:
- `config/config.yaml`: Configuration settings
- `workflow/`
  - `Snakefile`: Workflow definition
  - `envs/*.yaml`: Conda environment specifications
  - `rules/`: Workflow rules
  - `scripts/`: Custom Python scripts
  - `schemas/`: Configuration validation schemas
