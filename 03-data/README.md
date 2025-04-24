# Reference Data

This directory contains reference data and databases used in the analysis pipelines.

## Contents

- `all_50_bacterial_decoy.fa.tar.gz`: Bacterial decoy sequences used in mapping workflows
- `viruses_decoy.fasta.tar.gz`: Viral sequences found after taxonomy profile + bacterial decoy sequences
- `viruses.fasta.tar.gz`: Viral sequences found after taxonomy profile

## Database Organization

The reference data is organized to support various analysis tools:
- Reference sequences for mapping-based analyses
- Decoy sequences for filtering and quality control

## Usage

These reference files are used by various workflows in the `02-scripts` directory. The paths to these references should be configured in the respective workflow configuration files (`config.yaml`).
