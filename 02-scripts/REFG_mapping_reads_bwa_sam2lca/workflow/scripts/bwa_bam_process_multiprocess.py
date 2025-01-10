import pysam
import argparse
import sys
from multiprocessing import Pool
from itertools import islice
import os

########################################################################################

tab = "WSATUGCYRKMBDHVNwsatugcyrkmbdhvn".maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")

########################################################################################

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

########################################################################################

def revcomp(seq):
    return seq.translate(tab)[::-1]

########################################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description="Convert a BWA SAMSE/SAMPE BAM file with XA tags to Bowtie2-style BAM with separate secondary alignments.")

    parser.add_argument(
        "--input_bam", 
        default=snakemake.input.bam,
        type=str, 
        help="Input BAM file"
    )

    parser.add_argument(
        "--output_bam", 
        default=snakemake.output.bam,
        type=str, 
        help="Output BAM file with secondary alignments marked"
    )

    parser.add_argument(
        "-t",
        "--threads", 
        type=int, 
        default=snakemake.threads, 
        help="Number of threads for parallel processing"
    )

    parser.add_argument(
        "--chunk_size",
        type=int,
        default=10000,
        help="Number of reads to process in each chunk"
    )

    return parser.parse_args()
    
########################################################################################

def process_read(read, header):
    """
    Process a single read to split XA tag into secondary alignments.
    
    Args:
        read: The AlignedSegment object
        header: The BAM header
        
    Returns:
        List of tuples containing the data needed to recreate alignments
    """
    alignments = [(
        read.query_name,
        read.query_sequence,
        read.query_qualities,
        read.reference_id,
        read.reference_start,
        read.cigarstring,
        read.mapping_quality,
        read.flag,
        read.tags
    )]

    # Check if the XA tag exists
    if not read.has_tag("XA"):
        return alignments

    xa_tag = read.get_tag("XA")
    
    # Parse the XA tag
    alt_alignments = xa_tag.split(";")[:-1]  # Remove the trailing empty string after last semicolon

    for alt in alt_alignments:
        try:
            # Parse alternative alignment (chr,pos,CIGAR,NM)
            chr_name, pos, cigar, nm = alt.split(",")
            chr_name = chr_name.lstrip('+-')
            
            # Get reference ID from chromosome name
            ref_id = header.get('SQ', {}).get(chr_name, {}).get('id')
            if ref_id is None:
                continue
            
            # Calculate flag
            flag = read.flag | 0x100  # Set secondary alignment flag
            if pos.startswith('-'):
                flag |= 0x10  # Set reverse complement flag

            # Check if the read is reverse complemented
            if read.is_reverse != (flag & 0x10 != 0):
                sequence = revcomp(read.query_sequence)
                qualities = read.query_qualities[::-1]
            else:
                sequence = read.query_sequence
                qualities = read.query_qualities

            # Create alignment data tuple
            alignment_data = (
                read.query_name,
                sequence,
                qualities,
                ref_id,
                abs(int(pos)) - 1,  # Convert to 0-based
                cigar,
                0,  # mapping quality for secondary alignments
                flag,
                [('NM', int(nm))]  # Include NM tag
            )
            
            alignments.append(alignment_data)
            
        except Exception as e:
            print(f"Warning: Could not process alternative alignment {alt}: {str(e)}", file=sys.stderr)
            continue

    return alignments
########################################################################################

def process_chunk(args):
    """
    Process a chunk of reads.
    
    Args:
        args: Tuple containing (chunk of reads data, input BAM filename, header dict)
        
    Returns:
        List of alignment data tuples
    """
    chunk_data, input_bam, header = args
    all_alignments = []
    
    # Open input BAM to get header information
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        for read_data in chunk_data:
            # Recreate read from data
            read = pysam.AlignedSegment(bam_in.header)
            (read.query_name, read.query_sequence, read.query_qualities, 
             read.reference_id, read.reference_start, read.cigarstring,
             read.mapping_quality, read.flag, tags) = read_data
            for tag, value in tags:
                read.set_tag(tag, value)
            
            # Process read
            alignments = process_read(read, header)
            all_alignments.extend(alignments)
    
    return all_alignments
########################################################################################

def chunks(iterable, size):
    """Yield successive chunks from iterable"""
    iterator = iter(iterable)
    return iter(lambda: list(islice(iterator, size)), [])
########################################################################################

def process_bam(input_bam, output_bam, threads, chunk_size):
    # Open input BAM
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    header_dict = bam_in.header.to_dict()
    
    # Create a dictionary mapping reference names to their IDs
    ref_to_id = {name: idx for idx, name in enumerate(bam_in.references)}
    header_dict['SQ'] = {sq['SN']: {'id': ref_to_id[sq['SN']]} for sq in header_dict['SQ']}
    
    # Read all alignments and convert to serializable format
    print("Reading and chunking input BAM file...", file=sys.stderr)
    read_data = []
    count = 0
    for read in bam_in.fetch(until_eof=True):
        read_data.append((
            read.query_name,
            read.query_sequence,
            read.query_qualities,
            read.reference_id,
            read.reference_start,
            read.cigarstring,
            read.mapping_quality,
            read.flag,
            read.get_tags()
        ))
        count += 1
        if count % 100000 == 0:
            print(f"Read {count:,} alignments", file=sys.stderr)
    
    # Close input BAM
    bam_in.close()
    
    # Process chunks in parallel
    print(f"Processing {count:,} reads in parallel using {threads} threads...", file=sys.stderr)
    chunk_args = [(chunk, input_bam, header_dict) for chunk in chunks(read_data, chunk_size)]
    
    with Pool(processes=threads) as pool:
        all_results = []
        for i, chunk_results in enumerate(pool.imap_unordered(process_chunk, chunk_args)):
            all_results.extend(chunk_results)
            if (i + 1) % 10 == 0:
                print(f"Processed {(i + 1) * chunk_size:,} reads", file=sys.stderr)
    
    # Write results to output BAM
    print("Writing results to output BAM file...", file=sys.stderr)
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_in.header, threads=threads) as bam_out:
            for result in all_results:
                # Create new alignment
                new_read = pysam.AlignedSegment(bam_in.header)
                (new_read.query_name, new_read.query_sequence, new_read.query_qualities,
                 new_read.reference_id, new_read.reference_start, new_read.cigarstring,
                 new_read.mapping_quality, new_read.flag, tags) = result
                
                # Set tags
                for tag, value in tags:
                    new_read.set_tag(tag, value)
                
                bam_out.write(new_read)
    
    print(f"Processed BAM file saved as '{output_bam}'.", file=sys.stderr)
########################################################################################

def main(args):
    # Redirect stderr and stdout to log file if running in Snakemake
    if 'snakemake' in globals():
        sys.stderr = sys.stdout = open(snakemake.log[0], "w")
    
    # Run the conversion
    process_bam(args.input_bam, args.output_bam, args.threads, args.chunk_size)

########################################################################################

if __name__ == "__main__":
    args = parse_arguments()
    main(args)