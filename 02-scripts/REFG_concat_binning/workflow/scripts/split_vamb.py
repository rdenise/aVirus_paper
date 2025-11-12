import sys
import argparse
from Bio import SeqIO

########################################################################################

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

########################################################################################

num2sample = {
    "S1C": "ASM001",
    "S2C": "ASM002",
    "S3C": "ASM003",
    "S4C": "BSM001",
    "S5C": "BSM002",
    "S6C": "HSM001",
    "S7C": "HSM002",
    "S8C": "HSM003",
    "S9C": "HSM004",
    "S10C": "ZSM005",
    "S11C": "ZSM025",
    "S12C": "ZSM027",
    "S13C": "ZSM028",
    "S14C": "ZSM031",
    "S15C": "ZSM101",
    "S16C": "ZSM102",
    "S17C": "ZSM103",
    "S18C": "ZSM214",
    "S19C": "ZSM216",
    "S20C": "ZSM219",
}

sample2num = {
    "ASM001": "S1C",
    "ASM002": "S2C",
    "ASM003": "S3C",
    "BSM001": "S4C",
    "BSM002": "S5C",
    "HSM001": "S6C",
    "HSM002": "S7C",
    "HSM003": "S8C",
    "HSM004": "S9C",
    "ZSM005": "S10C",
    "ZSM025": "S11C",
    "ZSM027": "S12C",
    "ZSM028": "S13C",
    "ZSM031": "S14C",
    "ZSM101": "S15C",
    "ZSM102": "S16C",
    "ZSM103": "S17C",
    "ZSM214": "S18C",
    "ZSM216": "S19C",
    "ZSM219": "S20C"
}

########################################################################################

def split_fasta_by_sample(input_fasta, output_fasta):
    # Create dict to store sequences by sample
    seqs = []
    sample = output_fasta.split("/")[-1].split("-")[0]
    
    # Read input fasta and sort sequences by sample
    for record in SeqIO.parse(input_fasta, "fasta"):
        sample_id = record.id.split("C")[0].split("S")[1]  # Extract sample ID from header
        sample_id = f"S{sample_id}C"
        sample_name = num2sample.get(sample_id, None)

        if sample_name == sample:
            seqs.append(record)

    # Write sequences to sample-specific files
    if seqs:
        SeqIO.write(seqs, output_fasta, "fasta")
    else:
        with open(output_fasta, "w") as out_handle:
            pass

########################################################################################

def parse_args():
    parser = argparse.ArgumentParser(description='Split FASTA file by sample IDs')

    parser.add_argument(
        "-f", "--fasta", 
        default=snakemake.input.contigs,
        help="Contigs fasta file"
    )

    parser.add_argument(
        "-o", "--output",
        default=snakemake.output.fasta,
        help="Output file"
    )

    return parser.parse_args()

########################################################################################

if __name__ == "__main__":
    args = parse_args()
    split_fasta_by_sample(args.fasta, args.output)
