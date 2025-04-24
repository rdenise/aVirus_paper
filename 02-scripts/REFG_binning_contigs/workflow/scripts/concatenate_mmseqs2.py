import sys
import os
import argparse
import polars as pl

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

    parser.add_argument("--taxonomy_tsv", 
        help="Path to input taxonomy files",
        default=snakemake.input.tsv,
    )
    parser.add_argument("--output", 
        help="Paths to input TSV file", 
        default=snakemake.output.taxonomy,
    )

    # parse the arguments
    args = parser.parse_args()

    return args


###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()

    # Read the input TSV file
    list_df = []

    print(args.taxonomy_tsv)

    mapping = {
        ";Ceratitis;Ceratitis;": ";Ceratitis;",
        ";Lasius;Lasius;": ";Lasius;",
        ";Drosophila;Drosophila;": ";Drosophila;"
    }

    for tsv in args.taxonomy_tsv:
        df = pl.read_csv(tsv, separator="\t",
            has_header=False,
            new_columns=["contigs", "taxid", "rank", "name","predictions"],
        )
        
        df = df.with_columns(
            pl.col("predictions").str.replace_all(
                r"[a-z\-]_", ""
                ).str.split(
                    ";environmental samples"
                    ).list.first().str.replace_many(mapping)
        )

        df = df.select([
            "contigs",
            "predictions",
        ])
        
        list_df.append(df)

    # Concatenate the dataframes
    df = pl.concat(list_df)
    df.write_csv(args.output, separator="\t")
