# Analysis Results

This directory contains the results from both empirical and simulated data analyses.

## Directory Structure

### empirical_data/
Results from real sample analysis, organized by analysis type:
- `binning/`: Results from metagenomic binning
    - `semibin_bins.assembled.tsv.tar.gz`: Binning results using SemiBin with all assembled contigs
    - `vamb_bins.assembled.tsv.tar.gz`: Binning results using VAMB with all assembled contigs  
    - `taxvamb_bins.assembled.tsv.tar.gz`: Binning results using MMseqs2 with all assembled contigs
    - `semibin_bins.genomad.tsv.tar.gz`: Binning results using SemiBin with geNomad detected contigs
    - `vamb_bins.genomad.tsv.tar.gz`: Binning results using VAMB with geNomad detected contigs
    - `taxvamb_bins.genomad.tsv.tar.gz`: Binning results using MMseqs2 with geNomad detected contigs
    - `checkv_quality_summary_bins/`: CheckV quality assessment results for bins
        - `bins_from_assembled_contigs/`: Quality assessment for bins from all assembled contigs
            - `semibin_assembled.checkv_quality_summary.tsv.tar.gz`: CheckV results for SemiBin bins (assembled contigs)
            - `vamb_assembled.checkv_quality_summary.tsv.tar.gz`: CheckV results for VAMB bins (assembled contigs)
            - `taxvamb_assembled.checkv_quality_summary.tsv.tar.gz`: CheckV results for TaxVAMB bins (assembled contigs)
        - `bins_from_genomad_viral_contigs/`: Quality assessment for bins from geNomad detected viral contigs
            - `semibin_from_genomad_viral.checkv_quality_summary.tsv.tar.gz`: CheckV results for SemiBin bins (geNomad viral contigs)
            - `vamb_from_genomad_viral.checkv_quality_summary.tsv.tar.gz`: CheckV results for VAMB bins (geNomad viral contigs)
            - `taxvamb_from_genomad_viral.checkv_quality_summary.tsv.tar.gz`: CheckV results for TaxVAMB bins (geNomad viral contigs)
- `genomad_annotation/`: Viral genome detection and annotation results
    - `genomad.default.summary.virus.tsv.tar.gz`: Results from geNomad annotation (default mode, summary file)
    - `genomad.relaxed.summary.virus.tsv.tar.gz`: Results from geNomad annotation (relaxed mode, summary file)
    - `checkv_quality_summary_default.tsv.tar.gz`: Results from CheckV annotation (from genomad default mode)
    - `checkv_quality_summary_relaxed.tsv.tar.gz`: Results from CheckV annotation (from genomad relaxed mode)
    - `pydamage.filtered_results.genomad.default.tsv.tar.gz`: Results from pydamage analysis (filtered results, default mode)
    - `pydamage.filtered_results.genomad.relaxed.tsv.tar.gz`: Results from pydamage analysis (filtered results, relaxed mode)
- `mapping2refs/`: Results from reference mapping analyses
    - `combined_depth_breadth_read_contigs.tsv.tar.gz`: Combined depth and breadth of coverage for reads and contigs
    - `vclust_aln.tsv.tar.gz`: Results from vclust alignment file
    - `vclust_ani.tsv.tar.gz`: Results from vclust ani file
- `taxonomy_profile/`: Taxonomic classification results
    - `all.taxpasta.centrifuge.genus_viruses.txt.tar.gz`: Taxonomic classification using Centrifuge for all viruses at the genus level maximum
    - `all.taxpasta.centrifuge.species_viruses.txt.tar.gz`: Taxonomic classification using Centrifuge for all viruses at the species level maximum
    - `all.taxpasta.centrifuge.txt.tar.gz`: Taxonomic classification using Centrifuge for all viruses and bacterial decoy genomes

### simulation/
Results from simulated data analysis:
- `simulated_reads_taxonomy_profiler_comparison/`: Benchmarking results comparing different taxonomic profilers
    - `taxpasta.all.full.genus.virus.tsv.tar.gz`: Taxonomic classification results for all viruses at the genus level for all taxonomic profilers and truth (number of reads simulated)
    - `taxpasta.all.full.species.virus.tsv.tar.gz`: Taxonomic classification results for all viruses at the species level for all taxonomic profilers and truth (number of reads simulated)

## Data Organization

The analysis results are structured to facilitate:
- Comparison between different tools and methods
- Evaluation of tool performance using simulated data
- Integration of multiple lines of evidence for viral detection
- Assessment of ancient DNA characteristics

## Usage

These results can be used to:
- Compare the performance of different taxonomic profilers
- Evaluate the effectiveness of viral detection methods
- Assess the quality and authenticity of ancient DNA
- Generate figures and visualizations using scripts in the `02-scripts/PLOT_*` directories
