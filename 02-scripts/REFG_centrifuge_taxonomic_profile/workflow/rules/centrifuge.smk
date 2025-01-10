##########################################################################
##########################################################################

rule centrifuge_merge_reads:
    input:
        sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
    output:
        report = os.path.join(
            OUTPUT_FOLDER,
            "centrifuge",
            "{sample}.centrifuge.report.txt",
        ),
        output = os.path.join(
            OUTPUT_FOLDER,
            "centrifuge",
            "{sample}.centrifuge.output.txt",
        ),
        summary = os.path.join(
            OUTPUT_FOLDER,
            "centrifuge",
            "{sample}.centrifuge.summary.txt",
        ),
    params:
        extra="",  # optional parameters
        centrifuge_db=CENTRIFUGE_DB,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "centrifuge",
            "centrifuge.{sample}.log",
        ),
    threads: 20  # Use at least two threads
    resources:
        mem=80,
        cpus=20,
    conda: 
        "../envs/centrifuge.yaml"
    shell:
        """
        centrifuge -x {params.centrifuge_db:q} \
        --threads {threads} \
        -S {output.output:q} \
        --report-file {output.summary:q} \
        -k 1 \
        -U {input.sample_fastqs[2]:q} \
        -1 {input.sample_fastqs[0]:q} \
        -2 {input.sample_fastqs[1]:q} &> {log:q}

        centrifuge-kreport -x {params.centrifuge_db:q} \
        {output.output:q} 1> {output.report:q} 2>> {log:q}
        """

##########################################################################
##########################################################################

# rule centrifuge:
#     input:
#         sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
#     output:
#         report = os.path.join(
#             OUTPUT_FOLDER,
#             "centrifuge",
#             "{sample}.centrifuge.report.txt",
#         ),
#         output = os.path.join(
#             OUTPUT_FOLDER,
#             "centrifuge",
#             "{sample}.centrifuge.output.txt",
#         ),
#         summary = os.path.join(
#             OUTPUT_FOLDER,
#             "centrifuge",
#             "{sample}.centrifuge.summary.txt",
#         ),
#     params:
#         extra="",  # optional parameters
#         centrifuge_db=CENTRIFUGE_DB,
#     log:
#         os.path.join(
#             OUTPUT_FOLDER,
#             "logs",
#             "centrifuge",
#             "centrifuge.{sample}.log",
#         ),
#     threads: 20  # Use at least two threads
#     resources:
#         mem=80,
#         cpus=20,
#     conda: 
#         "../envs/centrifuge.yaml"
#     shell:
#         """
#         centrifuge -x {params.centrifuge_db:q} \
#         --threads {threads} \
#         -S {output.output:q} \
#         --report-file {output.summary:q} \
#         -k 1 \
#         -1 {input.sample_fastqs[0]:q} \
#         -2 {input.sample_fastqs[1]:q} &> {log:q}

#         centrifuge-kreport -x {params.centrifuge_db:q} \
#         {output.output:q} 1> {output.report:q} 2>> {log:q}
#         """

##########################################################################
##########################################################################
