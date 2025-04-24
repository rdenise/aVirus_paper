##########################################################################
##########################################################################

rule mmseqs2:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{sample}-{software}.fasta"
        ),
    output:
        report = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{sample}-{software}",
            "{sample}-{software}_lca.tsv",
        ),
        tmp = directory(
            os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{sample}-{software}_tmp"
            )
        )
    params:
        mmseqs_db=MMSEQS_DB,
        prefix = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{sample}-{software}",
            "{sample}-{software}"
        ),
        folder = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{sample}-{software}"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "mmseqs",
            "{method}",
            "{sample}-{software}",
            "mmseqs.{sample}-{software}.log",
        ),
    threads: 4  # Use at least two threads
    resources:
        mem=80,
        cpus=4,
    conda: 
        "../envs/mmseqs2.yaml"
    shell:
        """
        mkdir -p {params.folder:q}

        mmseqs easy-taxonomy \
        {input.contigs:q} \
        {params.mmseqs_db:q} \
        {params.prefix:q} \
        {output.tmp:q} \
        --tax-lineage 1 \
        --threads {threads} &> {log:q}
        """

##########################################################################
##########################################################################

rule concat_mmseqs:
    input:
        tsv=lambda wildcards: [
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                f"{wildcards.method}",
                f"{wildcards.software}",
                "mmseqs",
                f"{sample}-{wildcards.software}",
                f"{sample}-{wildcards.software}_lca.tsv"
            ) for sample in SAMPLES if sample.startswith(f"{wildcards.site}")
        ]
    output:
        taxonomy=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{site}.mmseqs2.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "mmseqs",
            "{site}-{software}.log",
        )
    resources:
        mem=40,
        cpus=1,
    threads: 1
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/concatenate_mmseqs2.py"

##########################################################################
##########################################################################