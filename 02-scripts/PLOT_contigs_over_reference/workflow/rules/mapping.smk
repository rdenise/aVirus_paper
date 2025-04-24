##########################################################################
##########################################################################

rule plot_contigs_over_reference:
    input: 
        genomad_relaxed = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus_summary.tsv",
        ),
        genomad_default = os.path.join(
            OUTPUT_FOLDER,
            "genomad_default",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus_summary.tsv",
        ),
        aln = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.aln.tsv",
        ),
        ani = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        viruses = config["viruses_contigs"],
    output:
        done = os.path.join(
            OUTPUT_FOLDER,
            "plot_contigs_over_reference",
            "{sample}",
            "{software}",
            "{sample}.{software}.done",
        )
    params:
        outdir = os.path.join(
            OUTPUT_FOLDER,
            "plot_contigs_over_reference",
            "{sample}",
            "{software}",
        ),
        sample = "{sample}",
        software = "{software}",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "plot_contigs_over_reference",
            "logs",
            "{sample}.{software}.log",
        ),
    wildcard_constraints:
        sample="[A-Z0-9]+",
        software="[a-zA-Z]+",
    resources:
        mem=50,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/plot_contigs_on_reference.py"

##########################################################################
##########################################################################
