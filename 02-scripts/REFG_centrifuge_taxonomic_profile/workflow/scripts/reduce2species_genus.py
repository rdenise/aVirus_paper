
import argparse
import sys
import os
from Bio import SeqIO
import gzip
import pandas as pd

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def reduce_genus(row):
    index_genus = row['rank_lineage'].split(';').index('genus')
    return ';'.join(row['lineage'].split(';')[:index_genus+1])

###########################################################
###########################################################

def reduce_species(row):
    index_species = row['rank_lineage'].split(';').index('species')
    return ';'.join(row['lineage'].split(';')[:index_species+1])

###########################################################
###########################################################

def parse_arguments():
    # create the argument parser
    parser = argparse.ArgumentParser(
        description="Reduce the contigs to only be in a size phage possible"
    )
    parser.add_argument(
        "--input_file",
        help="the path to the fasta taxpasta file",
        default=snakemake.input.tsv,
    )
    parser.add_argument(
        "--output_species",
        help="the path to the output file",
        default=snakemake.output.species,
    )
    parser.add_argument(
        "--output_genus",
        help="the path to the output file",
        default=snakemake.output.genus,
    )

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################


def main(args):

    # Read the input file
    taxpasta = pd.read_csv(
        args.input_file,
        sep='\t'
    ).fillna('')

    # Reduce to species
    condition = (taxpasta['rank_lineage'].str.contains('species$', regex=True)) | (taxpasta['rank_lineage'].str.contains('species;'))
    taxpasta.loc[condition, 'lineage'] = taxpasta[condition].apply(reduce_species, axis=1)

    previous_lineage = ''
    previous_index = 0
    records = taxpasta.to_dict(orient='records')
    index2remove = []

    for i, record in enumerate(records):
        if record['lineage'] == previous_lineage and previous_index != 0:
            taxpasta.iloc[previous_index,6:] += taxpasta.iloc[i,6:]
            index2remove.append(i)
        else:
            previous_lineage = record['lineage']
            previous_index = i

    taxpasta = taxpasta.drop(index2remove).reset_index(drop=True)

    taxpasta = taxpasta.rename(lambda column_name: column_name.split('.')[0], axis=1)

    taxpasta.to_csv(
        args.output_species,
        sep="\t",
        index=False,
    )

    taxpasta[taxpasta['lineage'].str.contains('Viruses')].to_csv(
        args.output_species.replace('.txt', '_viruses.txt'),
        sep="\t",
        index=False,
    )

    # Reduce to genus
    condition = taxpasta['rank_lineage'].str.contains('genus')
    taxpasta.loc[condition, 'lineage'] = taxpasta[condition].apply(reduce_genus, axis=1)

    previous_lineage = ''
    previous_index = 0
    records = taxpasta.to_dict(orient='records')
    index2remove = []

    for i, record in enumerate(records):
        if record['lineage'] == previous_lineage and previous_index != 0:
            taxpasta.iloc[previous_index,6:] += taxpasta.iloc[i,6:]
            index2remove.append(i)
        else:
            previous_lineage = record['lineage']
            previous_index = i

    taxpasta = taxpasta.drop(index2remove).reset_index(drop=True)

    taxpasta.to_csv(
        args.output_genus,
        sep="\t",
        index=False,
    )

    taxpasta[taxpasta['lineage'].str.contains('Viruses')].to_csv(
        args.output_genus.replace('.txt', '_viruses.txt'),
        sep="\t",
        index=False,
    )

    return

###########################################################
###########################################################


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
