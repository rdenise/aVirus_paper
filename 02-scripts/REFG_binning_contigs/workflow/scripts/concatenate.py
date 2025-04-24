#!/usr/bin/env python

import sys
import os
import argparse
import gzip
import vamb

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Creates the input FASTA file for Vamb.
    Input should be one or more FASTA files, each from a sample-specific assembly.
    If keepnames is False, resulting FASTA can be binsplit with separator 'C'.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )

    parser.add_argument("--outpath", 
        help="Path to output FASTA file",
        default=snakemake.output.fasta,
    )
    parser.add_argument("--inpaths", 
        help="Paths to input FASTA file(s)", 
        default=snakemake.input.refs,
    )
    parser.add_argument(
        "-m",
        dest="minlength",
        metavar="",
        type=int,
        default=400,
        help="Discard sequences below this length [2000]",
    )
    parser.add_argument(
        "--keepnames", 
        action="store_true", 
        help="Do not rename sequences [False]",
        default=False,
    )
    parser.add_argument(
        "--nozip", 
        action="store_true", 
        help="Do not gzip output [False]",
        default=False,
    )

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()

    # Run the code. Compressing DNA is easy, this is not much bigger than level 9, but
    # many times faster
    filehandle = (
        open(args.outpath, "w")
        if args.nozip
        else gzip.open(args.outpath, "wt", compresslevel=1)
    )

    print(args.inpaths)

    vamb.vambtools.concatenate_fasta(
        filehandle, args.inpaths, minlength=args.minlength, rename=(not args.keepnames)
    )
    filehandle.close()