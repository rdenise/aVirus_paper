#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict
import gzip
import polars as pl
from Bio import SeqIO, bgzf

###########################################################
# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Creates the input FASTA file for Vamb.
    Split the abundance table and concatenated file to be able to
    process by vamb.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    parser.add_argument("--abundance", 
        help="Path to input abundance file",
        default=snakemake.input.abundance,
    )
    parser.add_argument("--contigs", 
        help="Paths to input FASTA file(s)", 
        default=snakemake.input.concatenated,
    )
    parser.add_argument("--fasta_out", 
        help="Path to output FASTA files",
        nargs="+",
        default=snakemake.output.contigs,
    )
    parser.add_argument("--split_tsv", 
        help="Path to output abundance files", 
        nargs="+",
        default=snakemake.output.abundances,
    ),
    parser.add_argument("--samples", 
        help="Path to output FASTA files", 
        nargs="+",
        default=snakemake.params.samples,
    ),

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()

    print("Reading the abundance table and the concatenated contigs")

    print("The arguments are:")
    print("Abundance table:", args.abundance)
    print("Contigs:", args.contigs)
    print("Output FASTA files:", args.fasta_out)
    print("Output abundance files:", args.split_tsv)
    print("Samples:", args.samples)

    # Load the abundance table
    abundance = pl.read_csv(args.abundance, separator="\t")
    abundance_columns = abundance.columns

    num_samples = len(abundance_columns) - 1

    # Create a dictionary with the sample name and the output file
    dict_samples = defaultdict(str)
    i = 1
    
    end_file = os.path.basename(args.split_tsv[0]).split("-")[1]
    dirname = os.path.dirname(args.split_tsv[0])

    for sample in args.samples:
        site = f"S{i}C"

        sample = os.path.basename(sample)[:6]

        print(f"Sample {site} will be written to {sample}-{end_file}, in the directory {dirname}")

        dict_samples[site] = os.path.join(dirname, f"{sample}-{end_file}")
        i += 1

    # Split the abundance table
    for i in range(num_samples):
        site = f"S{i+1}C"

        output_abundance = dict_samples[site]
        print(f"Writing the abundance table for {site} to {output_abundance}")

        abundance.filter(
            pl.col('contigname').str.starts_with(site)
        ).write_csv(output_abundance, separator="\t")



    # As we will append to the file, we need to empty it before in case it exists
    dict_open_fasta = {site: open(dict_samples[site].replace(".abundance.tsv", ".fasta"), "w") for site in dict_samples}
        
    # Split the contigs

    with gzip.open(args.contigs, "rt") as handle:
        contigs = SeqIO.parse(handle, "fasta")

        for seq in contigs:
            site = seq.id.split("C")[0]
            site += "C"

            seq.name = seq.description = ""

            SeqIO.write(seq, dict_open_fasta[site], "fasta")

        for site in dict_open_fasta:
            dict_open_fasta[site].close()