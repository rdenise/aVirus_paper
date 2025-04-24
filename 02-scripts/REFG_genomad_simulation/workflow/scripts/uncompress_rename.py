
import argparse
import sys
import os
from collections import defaultdict
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
        "--input_file",
        help="the path to the fasta input files",
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
    dict_seq = {}
    with open(args.output_file, "w") as outfile:
        file = args.input_file
        if file.endswith(".gz"):
            with gzip.open(file, "rt") as handle:
                parser = SeqIO.parse(handle, 'fasta')

                for seq in parser:
                    seq.id = seq.id.replace(":+","_F").replace(":-","_R").replace(":","_").replace("-","_")
                    seq.name = seq.description = ""

                    if seq.id in dict_seq:
                        dict_seq[seq.id] += 1
                        seq.id = seq.id + "_" + str(dict_seq[seq.id])
                    else:
                        dict_seq[seq.id] = 1

                    SeqIO.write(seq, outfile, 'fasta')
        else:
            parser = SeqIO.parse(file, 'fasta')

            for seq in parser:
                seq.id = seq.id.replace(":+","_F").replace(":-","_R").replace(":","_").replace("-","_")
                seq.name = seq.description = ""

                if seq.id in dict_seq:
                    dict_seq[seq.id] += 1
                    seq.id = seq.id + "_" + str(dict_seq[seq.id])
                else:
                    dict_seq[seq.id] = 1
                    seq.id = seq.id + "_1"

                SeqIO.write(seq, outfile, 'fasta')

    
###########################################################
###########################################################


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
