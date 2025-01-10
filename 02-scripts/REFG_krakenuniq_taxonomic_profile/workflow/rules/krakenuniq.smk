##########################################################################
##########################################################################

rule krakenuniq_preload:
    input:
        sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
    output:
        output = os.path.join(
            OUTPUT_FOLDER,
            "krakenuniq",
            "{sample}.krakenuniq.output.txt",
        ),
        report = os.path.join(
            OUTPUT_FOLDER,
            "krakenuniq",
            "{sample}.krakenuniq.report.txt",
        ),
    params:
        extra="",  # optional parameters
        kraken_db=KRAKEN_DB,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "krakenuniq",
            "krakenuniq.{sample}.log",
        ),
    threads: 20  # Use at least two threads
    resources:
        mem=30,
        cpus=20,
    conda: 
        "../envs/krakenuniq.yaml"
    shell:
        """
        krakenuniq --db {params.kraken_db:q} \
        --preload \
        --threads {threads} \
        --output {output.output:q} \
        --report-file  {output.report:q} \
        --paired \
        {input.sample_fastqs:q} &> {log:q}
        """

##########################################################################
##########################################################################

# rule krakenuniq:
#     input:
#         sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
#     output:
#         output = os.path.join(
#             OUTPUT_FOLDER,
#             "krakenuniq",
#             "{sample}.krakenuniq.output.txt",
#         ),
#         report = os.path.join(
#             OUTPUT_FOLDER,
#             "krakenuniq",
#             "{sample}.krakenuniq.report.txt",
#         ),
#     params:
#         extra="",  # optional parameters
#         kraken_db=KRAKEN_DB,
#     log:
#         os.path.join(
#             OUTPUT_FOLDER,
#             "logs",
#             "krakenuniq",
#             "krakenuniq.{sample}.log",
#         ),
#     threads: 20  # Use at least two threads
#     resources:
#         mem=50,
#         cpus=20,
#     conda: 
#         "../envs/krakenuniq.yaml"
#     shell:
#         """
#         krakenuniq --db {params.kraken_db:q} \
#         --threads {threads} \
#         --output {output.output:q} \
#         --report-file  {output.report:q} \
#         --paired \
#         {input.sample_fastqs:q} &> {log:q}
#         """

##########################################################################
##########################################################################

# rule bracken:
#     input:
#         report = os.path.join(
#             OUTPUT_FOLDER,
#             "kraken2",
#             "{sample}.kraken2.report.txt",
#         ),
#     output:
#         report = os.path.join(
#             OUTPUT_FOLDER,
#             "bracken",
#             "{sample}.bracken.report.txt",
#         ),
#         output = os.path.join(
#             OUTPUT_FOLDER,
#             "bracken",
#             "{sample}.bracken.table.tsv",
#         ),
#     params:
#         extra="",  # optional parameters
#         kraken_db=KRAKEN_DB,
#     log:
#         os.path.join(
#             OUTPUT_FOLDER,
#             "logs",
#             "bracken",
#             "bracken.{sample}.log",
#         ),
#     threads: 1  # Use at least two threads
#     resources:
#         mem=15,
#         cpus=1,
#     conda: 
#         "../envs/kraken2.yaml"
#     shell:
#         """
#         bracken -d {params.kraken_db:q} \
#         -i {input.report:q} \
#         -o {output.output:q} \
#         -w {output.report:q} \
#         -r 75 \
#         -l S \
#         -t 0 &> {log:q}
#         """

##########################################################################
##########################################################################
