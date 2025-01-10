##########################################################################
##########################################################################

rule fastp_pe:
    input:
        sample=lambda wildcards: sample2file[wildcards.sample]
    output:
        trimmed=[
            os.path.join(
                OUTPUT_FOLDER,
                "fastq_merged",
                "{sample}.R1.fastq",
            ),
            os.path.join(
                OUTPUT_FOLDER,
                "fastq_merged",
                "{sample}.R2.fastq",
            ),
        ],
        # Unpaired reads separately
        unpaired1=os.path.join(
            OUTPUT_FOLDER,
            "fastq_merged",
            "{sample}.U1.fastq",
        ),
        unpaired2=os.path.join(
            OUTPUT_FOLDER,
            "fastq_merged",
            "{sample}.U2.fastq",
        ),
        # or in a single file
#        unpaired="trimmed/pe/{sample}.singletons.fastq",
        merged=os.path.join(
            OUTPUT_FOLDER,
            "fastq_merged",
            "{sample}.merged.fastq",
        ),
        failed=os.path.join(
            OUTPUT_FOLDER,
            "fastq_merged",
            "{sample}.failed.fastq",
        ),
        html=os.path.join(
            OUTPUT_FOLDER,
            "fastq_report",
            "{sample}.html",
        ),
        json=os.path.join(
            OUTPUT_FOLDER,
            "fastq_report",
            "{sample}.json",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fastp",
            "fastp.{sample}.log",
        ),
    params:
        extra="--merge -A -G -Q -L --overlap_len_require 11"
    threads: 2
    resources:
        mem=25,
        cpus=2
    wrapper:
        "v3.3.6/bio/fastp"

##########################################################################
##########################################################################
