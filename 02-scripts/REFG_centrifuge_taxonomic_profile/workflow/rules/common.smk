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
from collections import defaultdict

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
            "taxpasta",
            "all.taxpasta.{software}.species.txt",
        ),
        software=['centrifuge'],
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

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
FASTQ_FOLDER = config["metagenomes"]["reads_folder"]

sample2file = defaultdict(list)

separator = config["metagenomes"]["reads_identifier"]

# Real data
for fastqfile in glob.glob(os.path.join(FASTQ_FOLDER, "*.[Rm][120e]*")):
    sample = os.path.basename(fastqfile).split(separator)[0]
    sample2file[sample].append(fastqfile)

sample2file = {key: sorted(value) for key, value in sample2file.items() if key.startswith(("ASM", "BSM", "HSM", "ZSM"))}

# Simulated data
# for fastqfile in glob.glob(os.path.join(FASTQ_FOLDER, "*_[120]*")):
#     sample = os.path.basename(fastqfile).split(separator)[0]
#     sample2file[sample].append(fastqfile)

sample2file = {key: sorted(value) for key, value in sample2file.items()}

SAMPLE_NAMES = list(sample2file.keys())
    
CENTRIFUGE_DB = config["centrifuge_db"]

TAXONOMY_DB = config["taxonomy_db"]
