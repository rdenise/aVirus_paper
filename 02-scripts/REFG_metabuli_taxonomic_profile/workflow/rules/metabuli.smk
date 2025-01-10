##########################################################################
##########################################################################

rule metabuli:
    input:
        sample_fastqs=lambda wildcards: sample2file[wildcards.sample],
    output:
        classification = os.path.join(
            OUTPUT_FOLDER,
            "metabuli",
            "{sample}.metabuli_classifications.tsv",
        ),
        report = os.path.join(
            OUTPUT_FOLDER,
            "metabuli",
            "{sample}.metabuli_report.tsv",
        ),
        html = os.path.join(
            OUTPUT_FOLDER,
            "metabuli",
            "{sample}.metabuli_krona.html",
        ),
    params:
        extra="",  # optional parameters
        metabuli_db=METABULI_DB,
        outdir = os.path.join(
            OUTPUT_FOLDER,
            "metabuli",
        ),
        prefix = "{sample}.metabuli",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "metabuli",
            "metabuli.{sample}.log",
        ),
    threads: 20  # Use at least two threads
    resources:
        mem=150,
        cpus=20,
    conda: 
        "../envs/metabuli.yaml"
    shell:
        """
        metabuli classify \
        --threads {threads} \
        --min-score 0.15 --min-sp-score 0.5 \
        {input.sample_fastqs[0]:q} \
        {input.sample_fastqs[1]:q} \
        {params.metabuli_db:q} \
        {params.outdir:q} \
        {params.prefix} &> {log:q}
        """

##########################################################################
##########################################################################
