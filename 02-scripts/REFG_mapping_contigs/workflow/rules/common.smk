##########################################################################
##########################################################################
##
##                                Library
##
##########################################################################
##########################################################################

import os, sys
import glob
from snakemake.utils import validate
from Bio import SeqIO

##########################################################################
##########################################################################
##
##                               Functions
##
##########################################################################
##########################################################################


def get_final_output(outdir):
    """
    Generate final output name
    """
    final_output = []

    final_output += expand(
        os.path.join(
            outdir,
            "minimap2",
            "{sample}.contigs.{software}.sorted.bam",
        ),
        sample=SAMPLE_NAMES,
        software=["megahit", "metaspades"]
    )

    final_output += expand(
        os.path.join(
            outdir,
            "vclust",
            "{sample}.contigs.{software}.clusters.tsv",
        ),
        sample=SAMPLE_NAMES,
        software=["megahit", "metaspades"]
    )

    final_output += expand(
        os.path.join(
            outdir,
            "vclust",
            "depth_breadth",
            "{sample}.{software}.vclust.ani.depth_breadth.tsv",
        ),
        sample=SAMPLE_NAMES,
        software=["megahit", "metaspades"],
    )

    return final_output


##########################################################################


def create_folder(mypath):
    """
    Created the folder that I need to store my result if it doesn"t exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return


##########################################################################
##########################################################################
##
##                                Variables
##
##########################################################################
##########################################################################

# Validation of the config.yaml file
validate(config, schema="../schemas/config.schema.yaml")

##########################################################################
##########################################################################
##
##                           Options
##
##########################################################################
##########################################################################

# Result folder
OUTPUT_FOLDER = config["output_folder"]
# Adding to config for report
config["__output_folder__"] = os.path.abspath(OUTPUT_FOLDER)

## Multiple reference genomes
REFS_FILE = config["reference_genome"]

METADATA = config["metadata"]

CONTIGS_FOLDER = config["metagenomes"]["assemble_contigs"]

(SAMPLE_NAMES,) = glob_wildcards(
    os.path.join(CONTIGS_FOLDER, "megahit", "{sample_files}" + "-megahit.fasta.gz")
)

# SAMPLE_NAMES = [i for i in SAMPLE_NAMES if i.startswith(('ZSM103'))]
SAMPLE_NAMES = [i for i in SAMPLE_NAMES if i.startswith(('ASM', 'BSM', 'HSM', 'ZSM'))]

SAM2LCA_DB=config["sam2lca_db"]
SEQID2TAXID=config["seqid2taxid"]
TAXONOMY=config["taxonomy"]

