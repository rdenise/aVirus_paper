#########################################################
# 1. Importing libraries
#########################################################

import polars as pl
import matplotlib.pyplot as plt
from typing import List
from pathlib import Path
import sys
import argparse
from Bio import SeqIO

#########################################################
# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")


#########################################################
# 2. Functions
#########################################################

def plot_contig_matches(aln_df: pl.DataFrame, outfig: str, list_viral_relaxed: List[str], 
    list_viral_default: List[str], reference_id: str = "MG711460", reference_len: int = 36636):
    """
    Plot contig matches against a reference sequence
    
    Parameters:
    -----------
    aln_df : polars.DataFrame
       filter polars dataframe with columns query, reference, rstart, rend
    outfig : str
        Path to save the output figure
    list_viral_relaxed : List[str]
        Name of the contigs detected as viral in genomand relaxed
    list_viral_default : List[str]
        Name of the contigs detected as viral in genomand default
    reference_id : str
        ID of the reference sequence to plot matches for (default: MG711460)
    reference_len : int
        Length of the reference sequence (default: 36636)
    """
    
    custom_style = {
        "text.color": "#131516",
        "svg.fonttype": "none",
        "pdf.fonttype": 42
    }

    plt.rcParams.update(custom_style)
    
    # Add validation
    if aln_df.is_empty():
        return

    # Create plot
    with plt.style.context(custom_style):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), height_ratios=[1, 1])
        
        # Add reference lines
        ax1.hlines(y=-0.9, xmin=0, xmax=reference_len, color='grey', linestyle='-', linewidth=5)
        ax2.hlines(y=-0.9, xmin=0, xmax=reference_len, color='grey', linestyle='-', linewidth=5)
        
        # Split data into default and relaxed lists
        df = aln_df.filter(
            pl.col("reference").eq(reference_id),
            pl.col("query").str.starts_with("NODE_")
        )
        
        # Process default list (top panel)
        unique_queries = sorted(df['query'].unique())
        y_positions = {query: i for i, query in enumerate(unique_queries)}
        
        for row in df.to_dicts():
            y_pos = y_positions[row['query']]

            # Process default list (top panel)
            if row['query'] in list_viral_default:
                ax1.hlines(y_pos, row['rstart'], row['rend'], linewidth=5, color='green')
            else:
                ax1.hlines(y_pos, row['rstart'], row['rend'], linewidth=5, color='red')

            # Process relaxed list (bottom panel)
            if row['query'] in list_viral_relaxed:
                ax2.hlines(y_pos, row['rstart'], row['rend'], linewidth=5, color='green')
            else:
                ax2.hlines(y_pos, row['rstart'], row['rend'], linewidth=5, color='red')
        
        # Customize plots
        for ax, unique_queries, title in [(ax1, unique_queries, 'Default option'), 
                                        (ax2, unique_queries, 'Relaxed option')]:
            ax.grid(False)
            ax.set_ylabel('Contigs')
            ax.set_yticks(range(len(unique_queries)))
            ax.set_yticklabels(unique_queries)
            ax.set_title(f'{title} matches on {reference_id}')
        
        ax2.set_xlabel('Position on reference')
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(outfig, bbox_inches='tight')

        plt.close()
    return 

#########################################################
#########################################################

def parse_arguments():
    parser = argparse.ArgumentParser(description='Plot contigs on reference sequence')
    parser.add_argument(
        '--aln', 
        default=snakemake.input.aln, 
        help='Path to vclust TSV file'
    )
    parser.add_argument(
        '--ani', 
        default=snakemake.input.ani,
        help='Path to ANI TSV file'
    )
    parser.add_argument(
        '--genomad-default', 
        default=snakemake.input.genomad_default,
        help='Path to geNomad default TSV file'
    )
    parser.add_argument(
        '--genomad-relaxed', 
        default=snakemake.input.genomad_relaxed,
        help='Path to geNomad relaxed TSV file'
    )
    parser.add_argument(
        '--outdir', 
        default=snakemake.params.outdir,
        help='Path to output directory'
    )
    parser.add_argument(
        '--sample', 
        default=snakemake.params.sample,
        help='Sample name'
    )
    parser.add_argument(
        '--software', 
        default=snakemake.params.software,
        help='Software name'
    )
    parser.add_argument(
        '--output', 
        default=snakemake.output.done, 
        help='Output file'
    )
    parser.add_argument(
        '--viruses', 
        default=snakemake.input.viruses, 
        help='Path to viruses fasta file'
    )

    args = parser.parse_args()
    return args

#########################################################
# 3. Main
#########################################################

if __name__ == "__main__":
    args = parse_arguments()
    aln_tsv = args.aln
    ani_tsv = args.ani
    genomad_tsv = args.genomad_default
    genomad_relaxed_tsv = args.genomad_relaxed

    sample = args.sample
    software = args.software

    all_viruses = [i.id for i in SeqIO.parse(args.viruses, "fasta")]

    # Read blast results
    df = pl.read_csv(aln_tsv, separator="\t")

    df_ani = pl.read_csv(ani_tsv, separator="\t")
    reference_dict = df_ani.filter(
        pl.col("query").str.starts_with("NODE_"),
        ~pl.col("reference").str.starts_with("NODE_"),
        pl.col("reference").is_in(all_viruses)
        ).select(
            ['reference', "rlen"]
        ).rows_by_key('reference', unique=True)

    list_viral_default = pl.read_csv(genomad_tsv, separator="\t")
    list_viral_default = list_viral_default['seq_name'].unique().to_list()
    list_viral_default = [i.split("|")[0] for i in list_viral_default]

    list_viral_relaxed = pl.read_csv(genomad_relaxed_tsv, separator="\t")
    list_viral_relaxed = list_viral_relaxed['seq_name'].unique().to_list()
    list_viral_relaxed = [i.split("|")[0] for i in list_viral_relaxed]

    # Create output directory
    outfig = Path(args.outdir) 

    outfig.mkdir(parents=True, exist_ok=True)

    for reference_id, reference_len in reference_dict.items():        
        # Filter for reference sequence
        df_tmp = df.filter(pl.col("reference").eq(reference_id))

        # Create output file
        outfile = outfig / f"{sample}_{software}_{reference_id}.pdf"
        print(f"Plotting contigs on reference {reference_id} to {outfile}")

        # Plot contigs on reference
        p = plot_contig_matches(
            aln_df = df_tmp, 
            outfig = outfile, 
            list_viral_relaxed = list_viral_relaxed, 
            list_viral_default = list_viral_default, 
            reference_id = reference_id, 
            reference_len = reference_len,
        )

    with open(args.output, "w") as f:
        f.write("done")

