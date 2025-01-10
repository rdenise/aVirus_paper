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
            "{folder}plot_barplot",
            "simulated_assessment_{size}.1x.pdf",
        ),
        size=SIZE,
        folder=FOLDER,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "{folder}plot_sensitivity_vs_precision",
            "{size}.sensitivity_vs_precision.1x.pdf",
        ),
        size=SIZE,
        folder=FOLDER,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "{folder}plot_sensitivity_vs_precision",
            "sensitivity_vs_precision.1x.pdf",
        ),
        size=SIZE,
        folder=FOLDER,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "{folder}plot_{metric}_newsize",
            "{deamination}_deamination.{metric}_length.1x.pdf",
        ),
        deamination=DEAMINATION,
        folder=FOLDER,
        metric=["sensitivity", "precision"],
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

SIZE = ['Short', 'Medium', 'Long']
DEAMINATION = ['None', 'Light', 'Heavy']

PROFILER = ['krakenuniq', 'kraken2', 'centrifuge', 'metabuli']
# FOLDER = ['only_virus/', '']
FOLDER = ['']