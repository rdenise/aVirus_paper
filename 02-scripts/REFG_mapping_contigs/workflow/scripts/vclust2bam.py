#!/usr/bin/env python3
# vclust2bam.py

####################################################################################################
#
# Convert alignments to BAM with ANI and lengths
#
# This script reads a table with alignments and a table with ANI and lengths, and generates a BAM
# file with the alignments that pass a minimum ANI threshold.
#
# The input alignment table should have the following columns:
# - query: query sequence ID
# - reference: reference sequence ID
# - qlen: query sequence length
# - rlen: reference sequence length
# - alnlen: alignment length
# - rstart: reference start position
#
# The input ANI and lengths table should have the following columns:
# - query: query sequence ID
# - reference: reference sequence ID
#
# The output BAM file will have the following columns:
# - query: query sequence ID
# - flag: SAM flag (256 if the query is the same as the previous line, 0 otherwise)
# - reference: reference sequence ID
# - rstart: reference start position
# - mapq: mapping quality (ANI * 100)
# - cigar: CIGAR string (alignment length)
# - tlen: template length (query length)
#
# Usage: vclust2bam.py -i <input> -a <ani> -o <output> [-m <min_ani>]
#
# Required arguments:
# -i, --input: input alignment table
# -a, --ani: ANI and lengths file
# -o, --output: output BAM file
#
# Optional arguments:
# -m, --min-ani: minimum ANI threshold (default: 0.0)
#
####################################################################################################

import argparse
import polars as pl
import pysam
import sys
import tempfile
from pathlib import Path
import gzip
from Bio import SeqIO
import shutil
import os


########################################################################################

os.environ['POLARS_MAX_THREADS'] = "1"  # Control Polars threading
os.environ['OMP_NUM_THREADS'] = "1"  # Control OpenMP threading

####################################################################################################
#
# Functions
#
####################################################################################################


tab = "WSATUGCYRKMBDHVNwsatugcyrkmbdhvn".maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")

########################################################################################

def revcomp(seq: str) -> str:
    """
    Reverse complement a DNA sequence

    Parameters
    ----------
    seq : str
        DNA sequence
    
    Returns
    -------
    str
        Reverse complement of the input sequence
    """

    return seq.translate(tab)[::-1]

########################################################################################

def parse_args():
    """
    Parse command-line arguments
    """ 

    parser = argparse.ArgumentParser(description='Convert alignments to BAM with ANI and lengths')
    parser.add_argument(
        '-i', '--input', 
        default=snakemake.input.aln,
        help='Input alignment table',
    )
    parser.add_argument(
        '-a', '--ani', 
        default=snakemake.input.ani,
        help='ANI and lengths file',
    )
    parser.add_argument(
        '-o', '--output', 
        default=snakemake.output.bam,
        help='Output BAM file',
    )
    parser.add_argument(
        '-m', '--min-ani', 
        type=float, 
        default=snakemake.params.min_ani,
        help='Min ANI threshold',
    )
    parser.add_argument(
        '-f', '--fasta',
        default=snakemake.input.fasta,
        help='Input FASTA file (gzipped or not)',
    )
    parser.add_argument(
        '--metrics',
        default=snakemake.params.metrics,
        help='Choose the metrics to use (ani, tani, gani)',
        choices=['ani', 'tani', 'gani'],
    )

    return parser.parse_args()


####################################################################################################

def process_alignments(aln_file, ani_file, output_bam, min_ani, fasta_gz, metrics):
    """
    Process alignments and generate BAM file

    Parameters
    ----------
    aln_file : str
        Input alignment table
    ani_file : str
        ANI and lengths file
    output_bam : str
        Output BAM file
    min_ani : float
        Minimum ANI threshold
    fasta_gz : str
        Input FASTA file (gzipped)

    """

    # Read tables
    print(f"Reading {ani_file}")
    ani_df = pl.read_csv(ani_file, separator='\t')
    print(f"Reading {aln_file}")
    aln_df = pl.read_parquet(aln_file)
    
    if aln_df.shape[0] == 0:
        print("No alignments found")
        with open(output_bam, 'w') as f:
            f.write("@HD\tVN:1.0\tSO:coordinate\n")
        return

    # Read metadata
    reference_list = aln_df['reference'].unique().to_list()

    # Reduce ANI table to references in alignment table 
    ani_df = ani_df.filter(pl.col('reference').is_in(reference_list))

    # Join and process
    print("Joining and filtering")
    df = (aln_df.join(
        ani_df,
        left_on=['query', 'reference'],
        right_on=['query', 'reference']
    ).filter(
        pl.col(metrics) >= min_ani, 
        (pl.col('qcov').gt(0.85) | pl.col('rcov').gt(0.85)),
        pl.col('alnlen').gt(100),
        pl.col('reference').is_in(reference_list),
        pl.col('query').str.starts_with('NODE_')
    ).sort(['query', 'qcov', metrics, 'alnlen'], descending=[False, True, True, True])
    )

    stem_name = Path(output_bam).stem
    dirname = Path(output_bam).parent
    # Create temp directory and files
    print("Creating temporary files")
    with tempfile.TemporaryDirectory(prefix=f'tmp_{stem_name}', dir=dirname, ignore_cleanup_errors=True) as tmp_dir:
        tmp_dir = Path(tmp_dir)
        temp_sam = tmp_dir / f"{stem_name}.tmp.sam"
        temp_bam = tmp_dir / f"{stem_name}.tmp.bam"
        
        # Generate SAM header with reference sequences
        header = "@HD\tVN:1.0\tSO:coordinate\n"
        ref_lengths = ani_df.select(['reference', 'rlen']).unique()

        print(f"Writing header with {len(ref_lengths)} references")
        for row in ref_lengths.iter_rows(named=True):
            header += f"@SQ\tSN:{row['reference']}\tLN:{row['rlen']}\n"
        
        # Load sequences
        # fasta = Path(fasta_gz)
        # tmp_fasta = tmp_dir / fasta.stem

        print(f"Reading {fasta_gz}")
        # if fasta_gz.endswith('.gz'):
        #     with gzip.open(fasta_gz, 'rb') as f_in:
        #         with open(tmp_fasta, 'wb') as f_out:
        #             shutil.copyfileobj(f_in, f_out)

        #     sequences = SeqIO.index(str(tmp_fasta), 'fasta')
        # else:
        #     sequences = SeqIO.index(fasta_gz, 'fasta')
        sequences = SeqIO.index(fasta_gz, 'fasta')
        
        print(f"Writing {temp_sam}")
        with open(temp_sam, 'w') as sam_file:
            sam_file.write(header)
            current_query = None
            
            for row in df.iter_rows(named=True):
                flag = 256 if current_query == row['query'] else 0
                if flag == 0:
                    current_query = row['query']
                
                mapq = int(row['ani'] * 100)
                
                # Calculate soft clips
                left_clip = row['qstart'] - 1  # position is 1-based
                right_clip = row['qlen'] - (row['qstart'] + row['alnlen'] - 1)
                
                # Build CIGAR string from alignments
                cigar = ""
                current_op = None
                count = 0
                
                query = row['query_aln']
                ref = row['ref_aln']
                
                for q, r in zip(query, ref):
                    if q == '-':  # Deletion
                        if current_op == 'D':
                            count += 1
                        else:
                            if current_op:
                                cigar += f"{count}{current_op}"
                            current_op = 'D'
                            count = 1
                    elif r == '-':  # Insertion
                        if current_op == 'I':
                            count += 1
                        else:
                            if current_op:
                                cigar += f"{count}{current_op}"
                            current_op = 'I'
                            count = 1
                    else:  # Match/Mismatch
                        if current_op == 'M':
                            count += 1
                        else:
                            if current_op:
                                cigar += f"{count}{current_op}"
                            current_op = 'M'
                            count = 1
                
                if current_op:
                    cigar += f"{count}{current_op}"
                
                # Add soft clips if present
                if row['orientation'] == -1:
                    left_clip, right_clip = right_clip, left_clip
                    flag += 16
                
                if left_clip > 0:
                    cigar = f"{left_clip}S" + cigar
                if right_clip > 0:
                    cigar += f"{right_clip}S"
                
                tlen = row['qlen']
                
                # Get sequence and create dummy quality scores
                seq = str(sequences[row['query']].seq)

                if row['orientation'] == -1:
                    seq = revcomp(seq)

                qual = 'I' * len(seq)  # Using default Illumina quality score
                # qual = '*'  # Using a not applicable quality score
                
                sam_line = (f"{row['query']}\t{flag}\t{row['reference']}\t"
                           f"{row['rstart']}\t{mapq}\t{cigar}\t*\t0\t"
                           f"{tlen}\t{seq}\t{qual}\n")
                sam_file.write(sam_line)

        # Convert to sorted BAM
        print(f"Converting to sorted BAM {output_bam}")
        pysam.view("-b", "-o", str(temp_bam), str(temp_sam), catch_stdout=False)
        pysam.sort("-o", output_bam, str(temp_bam))
        pysam.index(output_bam)

    return

####################################################################################################
#
# Main
#
####################################################################################################

def main():
    args = parse_args()
    process_alignments(args.input, args.ani, args.output, args.min_ani, snakemake.input.fasta, args.metrics)

    stem_name = Path(args.output).stem
    dirname = Path(args.output).parent

    for tmp_dir in dirname.glob(f'tmp_{stem_name}*'):
        shutil.rmtree(tmp_dir)

####################################################################################################

if __name__ == "__main__":
    # Put error and out into the log file
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    main()