##########################################################################
##########################################################################

rule kraken2:
    input:
        sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
    output:
        report = os.path.join(
            OUTPUT_FOLDER,
            "kraken2",
            "{sample}.kraken2.report.txt",
        ),
        output = os.path.join(
            OUTPUT_FOLDER,
            "kraken2",
            "{sample}.kraken2.output.txt",
        ),
    params:
        extra="",  # optional parameters
        kraken_db=KRAKEN_DB,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "kraken2",
            "kraken2.{sample}.log",
        ),
    threads: 20  # Use at least two threads
    resources:
        mem=150,
        cpus=20,
    conda: 
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {params.kraken_db:q} \
        --threads {threads} \
        --output {output.output:q} \
        --report {output.report:q} \
        --confidence 0.15 \
        --use-names \
        --gzip-compressed \
        --paired \
        {input.sample_fastqs:q} &> {log:q}
        """

##########################################################################
##########################################################################

