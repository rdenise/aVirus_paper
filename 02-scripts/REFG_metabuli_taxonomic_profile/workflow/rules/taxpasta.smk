##########################################################################
##########################################################################

rule taxpasta:
    input:
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "{{software}}",
                "{sample}.{{software}}_report.tsv",
            ),
            sample = SAMPLE_NAMES,
        )
    output:
        tsv = os.path.join(
            OUTPUT_FOLDER,
            "taxpasta",
            "all.taxpasta.{software}.txt",
        ),
    params:
        extra="",  # optional parameters
        taxonomy=os.path.join(
            METABULI_DB,
            "taxonomy",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "taxpasta",
            "taxpasta.{software}.log",
        ),
    threads: 1  # Use at least two threads
    resources:
        mem=15,
        cpus=1,
    conda: 
        "../envs/taxpasta.yaml"
    shell:
        """
        taxpasta merge --profiler kraken2 \
        --output {output.tsv:q} \
        --output-format TSV \
        --wide \
        --taxonomy {params.taxonomy:q} \
        --add-name \
        --add-rank \
        --add-lineage \
        --add-id-lineage \
        --add-rank-lineage \
        {input:q} &> {log:q}
        """

##########################################################################
##########################################################################
