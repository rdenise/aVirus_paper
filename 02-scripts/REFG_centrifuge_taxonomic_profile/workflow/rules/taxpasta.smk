##########################################################################
##########################################################################

rule taxpasta:
    input:
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "{{software}}",
                "{sample}.{{software}}.report.txt",
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
        taxonomy=TAXONOMY_DB
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "taxpasta",
            "taxpasta.{software}.log",
        ),
    wildcard_constraints:
            software="[^.]+"
    threads: 1  # Use at least two threads
    resources:
        mem=15,
        cpus=1,
    conda: 
        "../envs/taxpasta.yaml"
    shell:
        """
        taxpasta merge --profiler centrifuge \
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

rule reduce2species_genus:
    input: 
        tsv = os.path.join(
                OUTPUT_FOLDER,
                "taxpasta",
                "all.taxpasta.{software}.txt",
            ),
    output: 
        species = os.path.join(
                OUTPUT_FOLDER,
                "taxpasta",
                "all.taxpasta.{software}.species.txt",
            ),
        genus = os.path.join(
                OUTPUT_FOLDER,
                "taxpasta",
                "all.taxpasta.{software}.genus.txt",
            ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "taxpasta",
            "taxpasta.reduce.{software}.log",
        ),
    resources:
        mem=40,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/reduce2species_genus.py" 

##########################################################################
##########################################################################
