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
            "sam2lca",
            "{sample}.{reference}.sam2lca.csv",
        ),
        sample=SAMPLE_NAMES,
        reference=REF_LIST,
    )

    final_output += expand(
        os.path.join(
            outdir,
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth.tsv",
        ),
        sample=SAMPLE_NAMES,
        reference=REF_LIST,
        rank=["genus", "species", "higher", "unclassified", "total"],
    )

    final_output += expand(
        os.path.join(
            outdir,
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth_cov10.tsv",
        ),
        sample=SAMPLE_NAMES,
        reference=REF_LIST,
        rank=["genus", "species", "higher", "total"],
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


REFS_FOLDER = Path(config["reference_genome"]).parent
REF_LIST = Path(config["reference_genome"]).stem

# path to contigs sheet (TSV format, columns: contig_name, path_contig)
FASTQ_FOLDER = config["metagenomes"]["reads_folder"]

sample2file = defaultdict(list)

separator = config["metagenomes"]["reads_identifier"]

for fastqfile in glob.glob(os.path.join(FASTQ_FOLDER, "*.[Rm][120e]*")):
    sample = os.path.basename(fastqfile).split(separator)[0]
    sample2file[sample].append(fastqfile)

sample2file = {key: sorted(value) for key, value in sample2file.items() if key.startswith(("ASM", "BSM", "HSM", "ZSM"))}

SAMPLE_NAMES = [i for i in list(sample2file.keys())]

SAM2LCA_DB=config["sam2lca_db"]
SEQID2TAXID=config["seqid2taxid"]
TAXONOMY=config["taxonomy"]