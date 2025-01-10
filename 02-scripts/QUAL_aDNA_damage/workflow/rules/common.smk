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

    # bam
    final_output += expand(
        os.path.join(
            outdir,
            "damage_profiling",
            "{sample}",
            "DamageProfiler.log",
        ),
        sample=SAMPLE_NAMES,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "pydamage",
            "{sample}",
            "pydamage_results.csv",
        ),
        sample=SAMPLE_NAMES,
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

REFS_FOLDER = config["reference_folder"]

REFS_FOLDER_VIRUS = config["ref_folder"]
# path to contigs sheet (TSV format, columns: contig_name, path_contig)
BAM_FOLDER = config["bam_folder"]

(SAMPLE_NAMES, SAMPLE_NAMES2,) = glob_wildcards(
        os.path.join(REFS_FOLDER, "{refs_folder}", "{sample}" + "refs.fasta")
    )

FASTQ_FOLDER = config["metagenomes"]["reads_folder"]

sample2file = defaultdict(list)

separator = config["metagenomes"]["reads_identifier"]

for fastqfile in glob.glob(os.path.join(FASTQ_FOLDER, "*_[12]*")):
    sample = os.path.basename(fastqfile).split(separator)[0]
    sample2file[sample].append(fastqfile)

sample2file = {key: sorted(value) for key, value in sample2file.items() if key.startswith(("ASM", "BSM", "HSM", "ZSM"))}
