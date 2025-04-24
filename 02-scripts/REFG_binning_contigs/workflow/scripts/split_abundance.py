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
    parser.add_argument("--fasta", 
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

def split_fasta(index_fasta, list_names, output_fasta):
    """
    Split a fasta file into a new fasta file with only the sequences
    in list_names.
    """

    with gzip.open(output_fasta, "wt") as handle:
        for id_contig in list_names:
            record = index_fasta[id_contig]
            record.name = record.description = ""
            SeqIO.write(record, handle, "fasta")

    return

###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()

    print("Reading the abundance table and the concatenated contigs")

    # Load the abundance table
    abundance = pl.read_csv(args.abundance, separator="\t")
    print(args.abundance)
    abundance_columns = abundance.columns

    abundance = abundance.with_columns(
        sites=pl.col("contigname").str.split("C").list[0].str.replace("S", "").cast(pl.Int32)
    )

    tmp_compress = args.contigs.replace("gz", "bgz")
    print(args.contigs)
    # Load the concatenated contigs
    with gzip.open(args.contigs, "rt") as handle:
        contigs = SeqIO.parse(handle, "fasta")
        
        with bgzf.BgzfWriter(tmp_compress, "wb") as outgz:
            SeqIO.write(contigs, outgz, "fasta")
    

    contigs = SeqIO.index(tmp_compress, "fasta")

    print("Splitting the abundance table and the concatenated contigs")
    
    dict_samples = defaultdict(list)
    i = 1

    print(args.samples)

    for sample in args.samples:
        site = os.path.basename(sample)[:3]
        dict_samples[site].append(i)
        i += 1

    print(dict_samples)

    print(args.fasta)

    for fasta_out in args.fasta:
        print(f"Processing {fasta_out}")
        site = os.path.basename(fasta_out).split(".")[0].split("-")[0]

        tmp_df = abundance.filter(
            pl.col("sites").is_in(dict_samples[site])
        )

        tmp_df = tmp_df.drop("sites")

        out_df = fasta_out.replace(".fna.gz", ".abundance.tsv")
        tmp_df.write_csv(out_df, separator="\t")

        list_contigs = tmp_df["contigname"].to_list()

        split_fasta(contigs, list_contigs, fasta_out)

    os.remove(tmp_compress)