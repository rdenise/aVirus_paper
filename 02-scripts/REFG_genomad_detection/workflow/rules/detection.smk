##########################################################################
##########################################################################

rule genomad:
    input:
        contigs=os.path.join(
            CONTIGS_FOLDER,
            "{software}",
            "{sample}-{software}.fasta.gz",
        ),
        db = GENOMAD_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        ),
    params:
        genomad_folder = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "genomad",
            "{sample}.genomad.{software}.log",
        ),
    resources:
        mem=10,
        cpus=10,
    threads: 10
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    conda: 
        "../envs/genomad.yaml"
    shell:
        """
        # genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
        genomad end-to-end --splits 8 --enable-score-calibration --restart --composition metagenome -t {threads} --cleanup --relaxed {input.contigs:q} {params.genomad_folder:q} {input.db:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule checkv:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        ),
        db = CHECKV_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "checkv",
            "{sample}-{software}",
            "quality_summary.tsv"
        ),
    params:
        checkv_folder = os.path.join(
            OUTPUT_FOLDER,
            "checkv",
            "{sample}-{software}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "{sample}.checkv.{software}.log",
        ),
    resources:
        mem=10,
        cpus=10,
    threads: 10
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    conda: 
        "../envs/checkv.yaml"
    shell:
        """
        # checkv end_to_end <input> <output> [options]
        checkv end_to_end {input.contigs} {params.checkv_folder} -t {threads} -d {input.db:q} --remove_tmp &> {log:q}
        """

##########################################################################
##########################################################################
