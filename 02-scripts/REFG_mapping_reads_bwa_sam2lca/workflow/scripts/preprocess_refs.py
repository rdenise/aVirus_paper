
import argparse
import sys
import os
from Bio import SeqIO
import gzip

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def parse_arguments():
    # create the argument parser
    parser = argparse.ArgumentParser(
        description="Reduce the contigs to only be in a size phage possible"
    )
    parser.add_argument(
        "--input_files",
        help="the path to the fasta input files",
        nargs="+",
        default=snakemake.input.refs,
    )
    parser.add_argument(
        "--output_file",
        help="the path to the output file",
        default=snakemake.output.fasta,
    )

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################


def main(args):
    with open(args.output_file, "w") as outfile:
        for file in args.input_files:
            if not os.path.isfile(file):
                raise FileNotFoundError(f"File {file} not found")
            
            if file.endswith(".gz"):
                with gzip.open(file, "rt") as handle:
                    parser = SeqIO.parse(handle, 'fasta')

                    for seq in parser:
                        SeqIO.write(seq, outfile, 'fasta')
            else:
                parser = SeqIO.parse(file, 'fasta')

                for seq in parser:
                    SeqIO.write(seq, outfile, 'fasta')

    
###########################################################
###########################################################


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
