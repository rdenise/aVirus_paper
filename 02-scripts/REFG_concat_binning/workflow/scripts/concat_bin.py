import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from typing import Dict, List, Tuple
from pathlib import Path
import sys
import tempfile
import gzip
import shutil
import os

########################################################################################

os.environ['POLARS_MAX_THREADS'] = "1"  # Control Polars threading
os.environ['OMP_NUM_THREADS'] = "1"  # Control OpenMP threading


# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

########################################################################################

def parse_args():
    """
    Parse command-line arguments
    """

    parser = argparse.ArgumentParser(description='Concat bins into a single sequence and save to fasta file')
    parser.add_argument(
        "-o", "--output",
        default=snakemake.output.fasta,
        help="Output file"
    )
    parser.add_argument(
        "-b", "--bins", 
        default=snakemake.params.clusters,
        help="Bins TSV file"
    )
    parser.add_argument(
        "-f", "--fasta", 
        default=snakemake.input.contigs,
        help="Contigs fasta file"
    )
    parser.add_argument(
        "-g", "--genomad", 
        default=snakemake.input.genomad,
        help="Genomad fasta file"
    )
    
    return parser.parse_args()

########################################################################################

def concat_bins(fasta_file: str, output_file: str, bins_file: str, genomad_file: str) -> None:
    """
    Concatenate bins into single sequences and write to fasta file

    Parameters
    ----------
    fasta_file : str
        Path to input fasta file
    output_file : str
        Path to output fasta file
    bins_file : str
        Path to bins TSV file
    genomad_file : str
        Path to genomad summary file
    """

    method = output_file.split("/")[-2]

    # Read bins
    if method == "semibin":
        bins_df = pl.read_csv(bins_file, separator="\t")
    else:
        bins_df = pl.read_csv(bins_file, separator="\t").rename({
            "clustername": "bin",
            "contigname": "contig"
        })

    # If genomad file exists, get viral contigs
    viral_contigs = set()
    if genomad_file:
        viral_contigs = SeqIO.parse(genomad_file, "fasta")
        viral_contigs = set(record.id.split("|")[0] for record in viral_contigs)

    # Group by bin
    bin_groups = bins_df.group_by("bin")

    # Create output fasta
    with open(output_file, "w") as out_handle:
        # Read input sequences
        sequences = SeqIO.index(fasta_file, "fasta")
        
        # Process each bin
        for bin_name, group in bin_groups:
            # Check if bin contains any viral contigs
            bin_contigs = set(group["contig"])
            if not viral_contigs or bin_contigs.intersection(viral_contigs):
                # Concatenate sequences
                combined_seq = ""
                for contig in bin_contigs:
                    if contig in sequences:
                        if combined_seq == "":
                            combined_seq = str(sequences[contig].seq)
                        else:
                            # Add a gap between sequences
                            combined_seq += "N" * 20
                            combined_seq += str(sequences[contig].seq)
                
                # Create new record
                record = SeqRecord(
                    Seq(combined_seq),
                    id=f"{method}_{bin_name[0]}",
                    description=""
                )
                
                # Write to file
                SeqIO.write(record, out_handle, "fasta")

########################################################################################

def main():
    """
    Main function
    """
    # Parse arguments
    args = parse_args()

    if not os.path.exists(args.bins):
        # Remove existing output file
        with open(args.output, "w") as f:
            pass
    else:
        # Concatenate bins
        if "assembled" in args.bins:
            genomad_file = args.genomad
        else:
            genomad_file = None
        concat_bins(args.fasta, args.output, args.bins, genomad_file)
    # Print completion message
    print(f"Concatenated bins saved to {args.output}")

########################################################################################

if __name__ == "__main__":
    # Run main function
    main()