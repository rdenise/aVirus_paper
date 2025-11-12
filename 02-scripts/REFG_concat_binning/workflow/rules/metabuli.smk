##########################################################################
##########################################################################

rule intall_taxconverter:
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "taxconverter.intall.done"
        )
    log: 
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "mmseqs2",
            "install_taxconverter.log"
        )
    conda: 
        "../envs/mmseqs2.yaml"
    shell: 
        """
        pip install git+https://github.com/RasmussenLab/taxconverter.git@main &> {log}
        touch {output}
        """

##########################################################################
##########################################################################

# rule split_abundance_metabuli:
#     input:
#         abundance=os.path.join(
#             OUTPUT_FOLDER,
#             "binning",
#             "{method}",
#             "{software}",
#             "strobealign",
#             "vamb",
#             "{software}.abundance.tsv"
#         ),
#         concatenated=os.path.join(
#             OUTPUT_FOLDER,
#             "binning",
#             "{method}",
#             "{software}",
#             "concatenation",
#             "{software}.concatenate.fasta.gz"
#         )
#     output:
#         abundances=expand(
#             os.path.join(
#                 OUTPUT_FOLDER,
#                 "binning",
#                 "{{method}}",
#                 "{{software}}",
#                 "strobealign",
#                 "vamb",
#                 "{sample}-{{software}}.abundance.tsv"
#             ),
#             sample=SAMPLES,
#         ),
#         contigs=expand(
#             os.path.join(
#                 OUTPUT_FOLDER,
#                 "binning",
#                 "{{method}}",
#                 "{{software}}",
#                 "strobealign",
#                 "vamb",
#                 "{sample}-{{software}}.fasta"
#             ),
#             sample=SAMPLES,
#         )
#     params:
#         samples=lambda wildcards: get_path_contigs(wildcards)
#     log:
#         os.path.join(
#             OUTPUT_FOLDER,
#             "logs",
#             "binning",
#             "{method}",
#             "strobealign",
#             "split",
#             "{software}.log",
#         )
#     resources:
#         mem=10,
#         cpus=1,
#     threads: 1
#     conda:
#         "../envs/python.yaml"
#     script: 
#         "../scripts/split_metabuli.py"

##########################################################################
##########################################################################

# checkpoint split_small_metabuli:
#     input:
#         fasta=os.path.join(
#             OUTPUT_FOLDER,
#             "binning",
#             "{method}",
#             "{software}",
#             "strobealign",
#             "vamb",
#             "{sample}-{software}.fasta"
#         ),
#     output:
#         contigs_dir=directory(
#             os.path.join(
#                 OUTPUT_FOLDER,
#                 "binning",
#                 "{method}",
#                 "{software}",
#                 "strobealign",
#                 "vamb",
#                 "split_contigs",
#                 "{sample}-{software}",
#             )
#         ),
#     log:
#         os.path.join(
#             OUTPUT_FOLDER,
#             "logs",
#             "binning",
#             "{method}",
#             "strobealign",
#             "split_small",
#             "{sample}-{software}.log",
#         )
#     resources:
#         mem=10,
#         cpus=1,
#     threads: 1
#     conda:
#         "../envs/python.yaml"
#     script: 
#         "../scripts/split_small.py"

##########################################################################
##########################################################################

rule metabuli:
    input:
        taxconverter = os.path.join(
            OUTPUT_FOLDER,
            "taxconverter.intall.done"
        ),
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "split_contigs",
            "{sample}-{software}",
            "{sample}-{software}.{sub}.fasta"
        ),
    output:
        classification = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "metabuli",
            "{sample}-{software}",
            "{sample}-{software}.{sub}.metabuli_classifications.tsv",
        ),
        report = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "metabuli",
            "{sample}-{software}",
            "{sample}-{software}.{sub}.metabuli_report.tsv",
        ),
        taxonomy = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "metabuli",
            "{sample}-{software}",
            "{sample}-{software}.{sub}.metabuli_taxonomy.tsv",
        ),
    params:
        metabuli_db=METABULI_DB,
        outdir = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "metabuli",
            "{sample}-{software}",
        ),
        prefix = "{sample}-{software}.{sub}.metabuli",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "metabuli",
            "{method}",
            "{sample}-{software}",
            "metabuli.{sample}-{software}.{sub}.log",
        ),
    threads: 4  # Use at least two threads
    resources:
        mem=40,
        cpus=4,
    conda: 
        "../envs/metabuli.yaml"
    shell:
        """
        metabuli classify \
        --seq-mode 3 \
        --threads {threads} \
        {input.contigs:q} \
        {params.metabuli_db:q} \
        {params.outdir:q} \
        {params.prefix} &> {log:q}

        taxconverter metabuli -c {output.classification} -r {output.report} -o {output.taxonomy} &>> {log}
        """

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

# def aggregate_input(wildcards):
#     # Get the checkpoint output
#     checkpoint_output = checkpoints.split_small_metabuli.get(**wildcards).output.contigs_dir

#     # Get the subsample names from the checkpoint output rule
#     subs = glob_wildcards(os.path.join(checkpoint_output, "{sample}-{software}.{sub}.fasta")).sub
    
#     # Return the list of files to aggregate from the metabuli rule
#     return [
#             os.path.join(
#                 OUTPUT_FOLDER,
#                 "binning",
#                 f"{wildcards.method}",
#                 f"{wildcards.software}",
#                 "metabuli",
#                 f"{wildcards.sample}-{wildcards.software}",
#                 f"{wildcards.sample}-{wildcards.software}.{sub}.metabuli_taxonomy.tsv",
#             )
#             for sub in subs
#         ]


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