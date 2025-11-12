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

    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "semibin",
    #         "{sample}-{software}.fasta",
    #     ),
    #     sample=SAMPLES,
    #     software=SOFTWARES,
    #     method=METHODS,
    # )   


    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "vamb",
    #         "{site}-{software}.fasta",
    #     ),
    #     software=SOFTWARES,
    #     method=METHODS,
    #     site=['ASM', 'BSM', 'HSM', 'ZSM'],
    # )

    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "taxvamb",
    #         "{site}-{software}.fasta",
    #     ),
    #     software=SOFTWARES,
    #     method=METHODS,
    #     site=['ASM', 'BSM'],
    # )

    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "taxvamb",
    #         "{sample}-{software}.fasta",
    #     ),
    #     software=SOFTWARES,
    #     method=METHODS,
    #     sample=[i for i in SAMPLES if i.startswith(('ZSM', 'HSM'))],
    # )

    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "{profiler}",
    #         "pydamage",
    #         "{sample}-{software}",
    #         "pydamage_filtered_results.csv",
    #     ),
    #     sample=SAMPLES,
    #     software=SOFTWARES,
    #     method=METHODS,
    #     profiler=["semibin", "vamb", "taxvamb"],
    # )  


    # final_output += expand(
    #     os.path.join(
    #         outdir,
    #         "binning",
    #         "{method}",
    #         "{software}",
    #         "concat_bins_fasta",
    #         "{profiler}",
    #         "checkv",
    #         "{sample}-{software}",
    #         "quality_summary.tsv"
    #     ),
    #     sample=SAMPLES,
    #     software=SOFTWARES,
    #     method=METHODS,
    #     profiler=["semibin", "vamb", "taxvamb"],
    # )  

    final_output += expand(
        os.path.join(
            outdir,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        sample=SAMPLES,
        software=SOFTWARES,
        method=METHODS,
        profiler=["semibin", "vamb", "taxvamb"],
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

for fastqfile in glob.glob(os.path.join(FASTQ_FOLDER, "*_[120]*")):
    sample = os.path.basename(fastqfile).split(separator)[0]
    sample2file[sample].append(fastqfile)

sample2file = {key: sorted(value) for key, value in sample2file.items()}

SAMPLES = list(sample2file.keys())

SAMPLES = [i for i in SAMPLES if i.startswith(('ASM', 'BSM', 'HSM', 'ZSM'))]

METABULI_DB = config["metabuli_db"]
MMSEQS_DB = config["mmseqs_db"]

METHODS = ["genomad", "assembled"]

SOFTWARES = ["megahit", "metaspades"]

REFS_FILE = config["reference_db"]