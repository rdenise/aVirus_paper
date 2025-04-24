import sys
import os
import argparse
from Bio import SeqIO

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

    parser.add_argument("--fasta", 
        help="Path to input abundance file",
        default=snakemake.input.fasta,
    )
    parser.add_argument("--outdir", 
        help="Paths to output directory", 
        default=snakemake.output.contigs_dir,
    )

    # parse the arguments
    args = parser.parse_args()

    return args

###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()

    threshold = 50000

    os.makedirs(args.outdir, exist_ok=True)

    parser = SeqIO.parse(args.fasta, "fasta")

    # Loop over records in chunks of 100
    seq_count = 1
    file_index = 1
    chunk = []

    sample_software = os.path.basename(args.fasta).split(".")[0]

    for seq in parser:
        seq_tmp = seq
        
        if len(seq_tmp.seq) <= threshold:
            seq_tmp.name = seq.description = ""

            chunk.append(seq_tmp)

            seq_count += 1

            if seq_count == 100:
                out_file = os.path.join(
                    args.outdir, 
                    f"{sample_software}.{file_index}.fasta"
                )
                
                SeqIO.write(chunk, out_file, "fasta")
                
                file_index += 1
                seq_count = 1
                chunk = []

    out_file = os.path.join(
        args.outdir, 
        f"{sample_software}.{file_index}.fasta"
    )
        

        