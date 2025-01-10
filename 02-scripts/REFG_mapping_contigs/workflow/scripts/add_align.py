import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from sequence_align.pairwise import hirschberg
from Bio.SeqRecord import SeqRecord
import argparse
import multiprocess as mp
from functools import partial
from typing import Dict, List, Tuple
from pathlib import Path
import sys
import tempfile
import gzip
import shutil
import os

########################################################################################

os.environ['POLARS_MAX_THREADS'] = "10"  # Control Polars threading
os.environ['OMP_NUM_THREADS'] = "1"  # Control OpenMP threading

########################################################################################

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

def align_sequences(seq1: str, seq2: str) -> Tuple[str, str]:
    """
    Align two sequences using the Hirschberg algorithm

    Parameters
    ----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence

    Returns
    -------
    Tuple[str, str]
        Aligned sequences
    """


    list_seq1 = list(seq1)
    list_seq2 = list(seq2)

    aligned_seq_a, aligned_seq_b = hirschberg(
        list_seq1,
        list_seq2,
        match_score=2,
        mismatch_score=-1,
        indel_score=-2   
    )

    return ''.join(aligned_seq_a), ''.join(aligned_seq_b)

########################################################################################

def parse_args():
    """
    Parse command-line arguments
    """

    parser = argparse.ArgumentParser(description='Add aligned sequences to alignment table')
    parser.add_argument(
        "-o", "--output",
        default=snakemake.output.parquet,
        help="Output file"
    )
    parser.add_argument(
        "-a", "--align", 
        default=snakemake.input.aln,
        help="Alignment TSV file"
    )
    parser.add_argument(
        "-m", "--metadata", 
        default=snakemake.input.metadata,
        help="Metadata TSV file"
    )
    parser.add_argument(
        "-r", "--ref", 
        default=snakemake.input.refs,
        help="Reference FASTA file"
    )
    parser.add_argument(
        "-q", "--query", 
        default=snakemake.input.contigs,
        help="Query FASTA file"
    )
    parser.add_argument(
        "-t", "--threads", 
        type=int, 
        default=snakemake.threads,
        help="Number of processes to use"
    )

    return parser.parse_args()

########################################################################################

def process_chunk(chunk_data: Tuple[int, pl.DataFrame], ref_seqs: Dict[str, SeqRecord], 
                 query_seqs: Dict[str, SeqRecord], tmp_dir: Path) -> Path:
    """Process chunk and save to temp parquet file"""
    chunk_id, chunk = chunk_data
    
    chunk_results = []
    
    for row in chunk.iter_rows(named=True):
        ref_id = row['reference']
        query_id = row['query']
        
        # Test if reference in + strand
        if row['rstart'] > row['rend']:
            rstart = row['rend']
            rend = row['rstart']
            ref_orient = -1 
        else:
            rstart = row['rstart']
            rend = row['rend']
            ref_orient = 1

        # Test if query in + strand
        if row['qstart'] > row['qend']:
            qstart = row['qend']
            qend = row['qstart']
            query_orient = -1
        else:
            qstart = row['qstart']
            qend = row['qend']
            query_orient = 1
        
        if ref_orient == query_orient:
            ref_seq = str(ref_seqs[ref_id].seq[rstart-1:rend])
            query_seq = str(query_seqs[query_id].seq[qstart-1:qend])
            query_orient = 1
        else:
            ref_seq = str(ref_seqs[ref_id].seq[rstart-1:rend])
            query_seq = revcomp(str(query_seqs[query_id].seq[qstart-1:qend]))
            query_orient = -1

        if qend - qstart != rend - rstart:
            aligned_ref, aligned_query = align_sequences(ref_seq, query_seq)
        else:
            aligned_ref, aligned_query = ref_seq, query_seq
            
        chunk_results.append({
            "query": query_id,
            "reference": ref_id,
            "pident": row['pident'],
            "alnlen": row['alnlen'],
            "qstart": qstart,
            "qend": qend,
            "rstart": rstart,
            "rend": rend,
            "orientation": query_orient,
            "nt_match": row['nt_match'],
            "nt_mismatch": row['nt_mismatch'],
            "query_aln": aligned_query,
            "ref_aln": aligned_ref
        })

    # Create and save chunk DataFrame
    chunk_df = pl.DataFrame(chunk_results)
    out_path = tmp_dir / f"chunk_{chunk_id}.parquet"
    chunk_df.write_parquet(out_path)
    
    return out_path

########################################################################################

def main(align_file: str, metadata_file: str, ref_fasta: str, query_fasta: str, n_processes: int, output: str):
    """
    Main function

    Parameters
    ----------
    align_file : str
        Alignment TSV file
    metadata_file : str
        Metadata TSV file
    ref_fasta : str
        Reference FASTA file
    query_fasta : str
        Query FASTA file
    n_processes : int
        Number of processes to use
    """

    # Load data
    print(f"Loading data: {align_file}", file=sys.stderr)
    df = pl.scan_csv(align_file, separator="\t").collect()

    print(f"Loading metadata: {metadata_file}", file=sys.stderr)
    metadata = pl.read_csv(metadata_file, separator="\t")
    metadata = metadata.rename({"Virus GENBANK accession": "reference"})
    
    print("Joining dataframes", file=sys.stderr)
    df = df.join(
        metadata,
        on="reference"
    ).filter(
        pl.col("query").str.starts_with("NODE_"),
        pl.col("Genus").str.contains('virus'),
    )

    # print(df.head(), file=sys.stderr)

    print(f"Loading reference sequences {ref_fasta}", file=sys.stderr)
    ref_seqs = SeqIO.index(ref_fasta, "fasta")
    stem_name = Path(align_file).stem
    dirname = Path(align_file).parent

    with tempfile.TemporaryDirectory(prefix=f'tmp_{stem_name}', dir=dirname, ignore_cleanup_errors=True) as tmp_dir:
        tmp_dir = Path(tmp_dir)
    #     with gzip.open(query_fasta, "rt") as f_in:
    #         tmp_fasta = tmp_dir / Path(query_fasta).stem

    #         print(f"Creating temporary directory: {tmp_dir} and unzipping {query_fasta} as {tmp_fasta}", file=sys.stderr)
    #         with open(tmp_fasta, "wt") as f_out:
    #             shutil.copyfileobj(f_in, f_out)
            
        # query_seqs = SeqIO.index(tmp_fasta, "fasta")
        query_seqs = SeqIO.index(query_fasta, "fasta")
        
        chunk_size = min(1000, df.shape[0] // n_processes + 1)
        chunks = [(i, chunk) for i, chunk in enumerate(df.iter_slices(n_rows=chunk_size))]
        
        print("Processing chunks", file=sys.stderr)
        try:
            with mp.get_context("spawn").Pool(n_processes) as pool:
                process_func = partial(process_chunk, ref_seqs=ref_seqs, query_seqs=query_seqs, tmp_dir=tmp_dir)
                chunk_files = pool.map(process_func, chunks, chunksize=1)
        except Exception as e:
            print(f"Multiprocessing error: {e}", file=sys.stderr)
            chunk_files = [process_func(chunk) for chunk in chunks]
            
        # Concatenate all chunk files
        print("Combining results", file=sys.stderr)
        # Print the first three lines of the files before reading them into polars
        if len(chunk_files) > 0:
            combined_df = pl.scan_parquet(chunk_files).collect()
            combined_df.write_parquet(output)
        else:
            with open(output, "w") as f:
                pl.DataFrame(
                    {
                        "query": [],
                        "reference": [],
                        "pident": [],
                        "alnlen": [],
                        "qstart": [],
                        "qend": [],
                        "rstart": [],
                        "rend": [],
                        "orientation": [],
                        "nt_match": [],
                        "nt_mismatch": [],
                        "query_aln": [],
                        "ref_aln": []
                    }
                ).write_parquet(f)
        
        # Clean up chunk files
        for f in chunk_files:
            f.unlink()

########################################################################################

if __name__ == "__main__":
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")
    args = parse_args()
    print("Starting", file=sys.stderr)
    main(args.align, args.metadata, args.ref, args.query, args.threads, args.output)