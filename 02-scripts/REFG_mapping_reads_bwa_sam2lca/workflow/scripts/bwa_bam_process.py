import pysam
import argparse
import sys

tab = "WSATUGCYRKMBDHVNwsatugcyrkmbdhvn".maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")

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

    return parser.parse_args()

########################################################################################


def process_read(read, bam_in):
    """
    Process a single read to split XA tag into secondary alignments.
    
    Args:
        read: The AlignedSegment object
        bam_in: The input BAM file object (needed for header information)
        
    Returns:
        List of AlignedSegment objects, including the primary and secondary alignments.
    """
    alignments = [read]

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
            strand, pos = pos[0], pos[1:]
            chr_name = chr_name.lstrip('+-')
            
            # Verify that the chromosome exists in the header
            if chr_name not in bam_in.references:
                print(f"Warning: Reference {chr_name} not found in BAM header, skipping alignment", file=sys.stderr)
                continue
            
            # Create a new AlignedSegment object with the same header
            sec_alignment = pysam.AlignedSegment(bam_in.header)
            
            # Copy basic attributes from the primary alignment
            sec_alignment.query_name = read.query_name
            sec_alignment.query_sequence = read.query_sequence
            sec_alignment.query_qualities = read.query_qualities
            
            # Set the alignment specific attributes
            sec_alignment.reference_id = bam_in.get_tid(chr_name)
            sec_alignment.reference_start = abs(int(pos)) - 1  # Convert to 0-based
            sec_alignment.cigarstring = cigar
            sec_alignment.mapping_quality = 0  # Secondary alignments typically have mapping quality 0
            sec_alignment.set_tag("NM", int(nm))

            # Set flags
            sec_alignment.is_secondary = True
            sec_alignment.is_reverse = strand == '-'

            # Set the sequence
            if sec_alignment.is_reverse != read.is_reverse:
                sec_alignment.query_sequence = revcomp(read.query_sequence)
                sec_alignment.query_qualities = read.query_qualities[::-1]

            
            # Copy paired-end flags if applicable
            if read.is_paired:
                sec_alignment.is_paired = True
                sec_alignment.is_read1 = read.is_read1
                sec_alignment.is_read2 = read.is_read2
            
            alignments.append(sec_alignment)
            
        except Exception as e:
            print(f"Warning: Could not process alternative alignment {alt}: {str(e)}", file=sys.stderr)
            continue

    return alignments

########################################################################################

def process_bam(input_bam, output_bam, threads):
    # Open input BAM
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    
    # Open output BAM with the same header as input
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header, threads=threads)

    # Process reads
    count = 0
    for read in bam_in.fetch(until_eof=True):
        alignments = process_read(read, bam_in)
        for alignment in alignments:
            bam_out.write(alignment)
        
        count += 1
        if count % 100000 == 0:
            print(f"Processed {count:,} reads", file=sys.stderr)

    # Close files
    bam_in.close()
    bam_out.close()
    print(f"Processed BAM file saved as '{output_bam}'.", file=sys.stderr)

########################################################################################


def main(args):
    # Redirect stderr and stdout to log file if running in Snakemake
    if 'snakemake' in globals():
        sys.stderr = sys.stdout = open(snakemake.log[0], "w")
    
    # Run the conversion
    process_bam(args.input_bam, args.output_bam, args.threads)

########################################################################################

if __name__ == "__main__":
    args = parse_arguments()
    main(args)