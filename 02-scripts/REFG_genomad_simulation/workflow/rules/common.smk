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

    # bam

    final_output += expand(
        os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "checkv_{software}",
            "{genome}",
            "quality_summary.tsv"
        ),
        genome=SAMPLE_NAMES,
        software=['genomad', 'genomad_default']
    )

    final_output += expand(
        os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "jeager",
            "split_contigs",
            "{genome}",
            "{genome}_default_phages_jaeger.fasta",
        ),
        genome=SAMPLE_NAMES,
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
    os.path.join(CONTIGS_FOLDER, "{sample_files}" + ".fna.gz")
)

GENOMAD_DB=config["genomad_db"]

CHECKV_DB=config["checkv_db"]




