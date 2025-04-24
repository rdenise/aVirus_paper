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
            "plot_contigs_over_reference",
            "{sample}",
            "{software}",
            "{sample}.{software}.done",
        ),
        sample=SAMPLES,
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

# Samples
(SAMPLES, SAMPLES2, ) = glob_wildcards(
    os.path.join(OUTPUT_FOLDER, "genomad", "{sample}-megahit", "{sample2}-megahit_summary.log")
)
