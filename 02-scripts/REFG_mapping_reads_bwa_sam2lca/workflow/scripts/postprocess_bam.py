
import argparse
import sys
import os
import multiprocessing
from functools import partial
import tempfile
from itertools import islice
import pysam

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description="Split BAM file by taxonomic tag.")
    parser.add_argument(
        "--input_bam", 
        type=str, 
        default=snakemake.input.bam, 
        help="Input BAM file")

    parser.add_argument(
        "--output_genus_bam", 
        type=str, 
        default=snakemake.output.genus, 
        help="Output BAM file for genus")

    parser.add_argument(
        "--output_species_bam", 
        type=str, 
        default=snakemake.output.species, 
        help="Output BAM file for species")

    parser.add_argument(
        "--output_other_bam", 
        type=str, 
        default=snakemake.output.other, 
        help="Output BAM file for other")

    parser.add_argument(
        "--output_unclassified_bam",
        type=str,
        default=snakemake.output.unclassified,
        help="Output BAM file for unclassified reads")

    parser.add_argument(
        "--threads",
        type=int,
        default=snakemake.threads,
        help="Number of threads to use")

    return parser.parse_args()

###########################################################
###########################################################

def split_bam_by_taxonomic_tag(input_bam, output_genus_bam, output_species_bam, output_other_bam, output_unclassified_bam):
    # Open the input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as bamfile:

        with pysam.AlignmentFile(output_genus_bam, "wb", header=bamfile.header) as genus_bam, \
            pysam.AlignmentFile(output_species_bam, "wb", header=bamfile.header) as species_bam, \
            pysam.AlignmentFile(output_unclassified_bam, "wb", header=bamfile.header) as unclassified_bam, \
            pysam.AlignmentFile(output_other_bam, "wb", header=bamfile.header) as other_bam:

            # Iterate through each read in the BAM file
            for read in bamfile:
                tags = dict(read.get_tags())

                # Check for XR tag and classify by genus or species
                if 'XR' in tags:
                    if tags['XR'].lower() == 'genus':
                        genus_bam.write(read)
                    elif tags['XR'].lower() == 'species':
                        species_bam.write(read)
                    elif tags['XR'].lower() == 'no rank' and tags['XT'] not in {1, 131567, 2646395}:
                        species_bam.write(read)
                    else:
                        other_bam.write(read)
                else:
                    unclassified_bam.write(read)

    return 

###########################################################
###########################################################

def get_total_reads(bam_file):
    """Get total number of reads in BAM file"""
    count = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for _ in bam:
            count += 1
    return count

###########################################################
###########################################################

def process_chunk(chunk_info, input_bam, temp_dir):
    """Process a chunk of reads specified by start and end indices"""
    start_idx, end_idx = chunk_info
    
    # Initialize temporary BAM files for this chunk
    temp_paths = {}
    temp_files = {}
    
    try:
        with pysam.AlignmentFile(input_bam, "rb") as input_bamfile:
            header = input_bamfile.header.copy()
            
            for category in ['genus', 'species', 'other', 'unclassified']:
                temp_path = os.path.join(temp_dir, f"temp_{category}_chunk_{start_idx}_{end_idx}.bam")
                temp_paths[category] = temp_path
                temp_files[category] = pysam.AlignmentFile(temp_path, "wb", header=header)
            
            # Skip to start position
            for _ in islice(input_bamfile, start_idx):
                pass
            
            # Process assigned chunk of reads
            for read in islice(input_bamfile, end_idx - start_idx):
                tags = dict(read.get_tags())
                
                if 'XR' in tags:
                    if tags['XR'].lower() == 'genus':
                        temp_files['genus'].write(read)
                    elif tags['XR'].lower() == 'species':
                        temp_files['species'].write(read)
                    elif tags['XR'].lower() == 'no rank' and tags['XT'] not in {1, 131567, 2646395}:
                        temp_files['species'].write(read)
                    else:
                        temp_files['other'].write(read)
                else:
                    temp_files['unclassified'].write(read)
    
    finally:
        # Close all temporary files
        for f in temp_files.values():
            f.close()
    
    return temp_paths

###########################################################
###########################################################

def split_bam_by_taxonomic_tag_parallel(input_bam, output_genus_bam, output_species_bam, output_other_bam, output_unclassified_bam, threads):
    # Get total number of reads
    total_reads = get_total_reads(input_bam)
    chunk_size = total_reads // threads + (1 if total_reads % threads else 0)
    
    # Create chunks based on read indices
    chunks = []
    for i in range(0, total_reads, chunk_size):
        end = min(i + chunk_size, total_reads)
        chunks.append((i, end))
    
    print(f"Processing {total_reads} reads in {len(chunks)} chunks using {threads} threads")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory(dir=os.path.dirname(output_genus_bam)) as temp_dir:
        # Initialize multiprocessing pool
        with multiprocessing.Pool(threads) as pool:
            # Process chunks in parallel
            process_func = partial(process_chunk, input_bam=input_bam, temp_dir=temp_dir)
            results = pool.map(process_func, chunks)
        
        # Collect all temporary file paths
        all_temp_files = {
            'genus': [],
            'species': [],
            'other': [],
            'unclassified': []
        }
        
        for result in results:
            for category in all_temp_files:
                if os.path.exists(result[category]):
                    all_temp_files[category].append(result[category])
        
        # Get header from input BAM
        with pysam.AlignmentFile(input_bam, "rb") as input_bamfile:
            header = input_bamfile.header.copy()
        
        # Merge results for each category
        if all_temp_files['genus']:
            pysam.merge("-f", output_genus_bam, *all_temp_files['genus'], "-@", str(threads))
        else:
            pysam.AlignmentFile(output_genus_bam, "wb", header=header).close()
            
        if all_temp_files['species']:
            pysam.merge("-f", output_species_bam, *all_temp_files['species'], "-@", str(threads))
        else:
            pysam.AlignmentFile(output_species_bam, "wb", header=header).close()
            
        if all_temp_files['unclassified']:
            pysam.merge("-f", output_unclassified_bam, *all_temp_files['unclassified'], "-@", str(threads))

        if all_temp_files['other']:
            pysam.merge("-f", output_other_bam, *all_temp_files['other'], "-@", str(threads))
        else:
            pysam.AlignmentFile(output_other_bam, "wb", header=header).close()

###########################################################
###########################################################

def main(args):
    # Check if number of threads is greater than 1
    if args.threads > 1:
        split_bam_by_taxonomic_tag_parallel(
            input_bam = args.input_bam, 
            output_genus_bam = args.output_genus_bam, 
            output_species_bam = args.output_species_bam, 
            output_other_bam = args.output_other_bam,
            output_unclassified_bam = args.output_unclassified_bam,
            threads = args.threads)
    else:
        split_bam_by_taxonomic_tag(
            input_bam = args.input_bam, 
            output_genus_bam = args.output_genus_bam, 
            output_species_bam = args.output_species_bam, 
            output_other_bam = args.output_other_bam, 
            output_unclassified_bam = args.output_unclassified_bam
        )


###########################################################
###########################################################

if __name__ == "__main__":
    args = parse_arguments()
    main(args)