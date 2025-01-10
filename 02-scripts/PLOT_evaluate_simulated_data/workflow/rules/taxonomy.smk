##########################################################################
##########################################################################

rule taxonomy_postprocess:
    input: 
        tsv = lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "{taxofile}.{profiler}.output.txt" if "metabuli" != wildcards.profiler else "{taxofile}.{profiler}_classifications.tsv",
        ),
        metadata = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/ICTV_database/20240209/ICTV_metadata.tsv',
    output:
        taxonomy = os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "{taxofile}.{profiler}.taxonomy.parquet",
        )
    params:
        taxo_profiler = "{profiler}",
        seq2name = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/krakenuniq_ICTV/seqid2taxid.map',
        seq2name_old = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/krakenuniq_ICTV/seqid2taxid.map.orig',
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "logs",
            "{taxofile}.{profiler}.log",
        ),
    wildcard_constraints:
        folder="[a-z0-9]+",
        taxofile="[a-z0-9-]+",
    resources:
        mem=70,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/taxonomy_reformat.py"

##########################################################################
##########################################################################

rule taxonomy_postprocess_only_virus:
    input: 
        tsv = lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "only_virus",
            "{profiler}",
            "{taxofile}.{profiler}.output.txt" if "metabuli" != wildcards.profiler else "{taxofile}.{profiler}_classifications.tsv"
        ),
        metadata = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/ICTV/ICTV_database/20240209/ICTV_metadata.tsv',
    output:
        taxonomy = os.path.join(
            OUTPUT_FOLDER,
            "only_virus",
            "{profiler}",
            "{taxofile}.{profiler}.taxonomy.parquet",
        )
    params:
        taxo_profiler = "{profiler}",
        seq2name = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/krakenuniq_db_base/virus_only/seqid2taxid.map',
        seq2name_old = '/mnt/archgen/microbiome_coprolite/aVirus/03-data/refdbs/krakenuniq_db_base/virus_only/seqid2taxid.map.orig',
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "only_virus",
            "{profiler}",
            "logs",
            "{taxofile}.{profiler}.log",
        ),
    resources:
        mem=70,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/taxonomy_reformat_only_virus.py"

##########################################################################
##########################################################################

rule barplot:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                OUTPUT_FOLDER,
                "{folder}",
                f"{wildcards.size.lower()}-lognormal-{deamination}-deamination-{coverage}.{{profiler}}.taxonomy.parquet",
            ) for coverage in [0, 1, 2, 3, 4, 5] for deamination in ["none", "light", "heavy"]
        ]
    output:
        tsv = os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "{profiler}.{size}.count_melt.parquet",
        ),
        pdf = os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "plot",
            "{profiler}_simulated_assessment_{size}.pdf",
        )
    params:
        taxo_profiler = "{profiler}",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "logs",
            "{size}.barplot.{profiler}.log",
        ),
    resources:
        mem=150,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/bar_plot_correct.py"

##########################################################################
##########################################################################

rule barplot_Alex:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                "{folder}",
                f"{profiler}",
                f"{wildcards.size.lower()}-lognormal-{deamination}-deamination-2.{profiler}.taxonomy.parquet", # 2 = 1x for phage
            ) for profiler in [
                "metabuli",
                "krakenuniq",
                "kraken2",
                "centrifuge",
            ] for deamination in ["none", "light", "heavy"]
        ]
    output:
        tsv = os.path.join(
            "{folder}",
            "plot_barplot",
            "{size}.count_melt.parquet",
        ),
        pdf = os.path.join(
            "{folder}",
            "plot_barplot",
            "simulated_assessment_{size}.1x.pdf",
        )
    log:
        os.path.join(
            "{folder}",
            "plot_barplot",
            "logs",
            "{size}.barplot.log",
        ),
    resources:
        mem=40,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/bar_plot_alex.py"

##########################################################################
##########################################################################

rule minimizers_vs_reads:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                OUTPUT_FOLDER,
                "{folder}",
                f"{wildcards.size.lower()}-lognormal-{deamination}-deamination-2.{{profiler}}.taxonomy.parquet", # 2 = 1x for phage
            ) for deamination in ["none", "light", "heavy"]
        ],
        inspect = lambda wildcards: os.path.join(
            "/mnt",
            "archgen",
            "microbiome_coprolite",
            "aVirus",
            "03-data",
            "refdbs",
            "ICTV",
            "kraken2_only_ICTV" if "only_virus" in wildcards.folder else "kraken2_ICTV",
            "kraken.inspect.only_ICTV.taxpasta.tsv" if "only_virus" in wildcards.folder else "kraken.inspect.ICTV.taxpasta.tsv",
        )
    output:
        tsv = os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "{profiler}.{size}.1x.correct_melt.parquet",
        ),
        pdf = os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "plot",
            "{profiler}_simulated_correct_annotation_{size}.1x.pdf",
        )
    params:
        taxo_profiler = "{profiler}",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "{folder}",
            "logs",
            "{size}.minimizer_vs_reads.{profiler}.log",
        ),
    resources:
        mem=40,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/minimizers_vs_reads.py"

##########################################################################
##########################################################################

rule sensitivity_vs_precision:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                "{folder}",
                f"{profiler}",
                f"{size}-lognormal-heavy-deamination-2.{profiler}.taxonomy.parquet", # 2 = 1x for phage
            ) for profiler in [
                "metabuli",
                "krakenuniq",
                "kraken2",
                "centrifuge",
            ] for size in ["short", "medium", "long"]
        ]
    output:
        pdf = os.path.join(
            "{folder}",
            "plot_sensitivity_vs_precision",
            "sensitivity_vs_precision.1x.pdf",
        ),
        png = os.path.join(
            "{folder}",
            "plot_sensitivity_vs_precision",
            "sensitivity_vs_precision.1x.png",
        )
    params:
        tsv = os.path.join(
            "{folder}",
            "plot_sensitivity_vs_precision",
            "sensitivity_precision.parquet",
        ),
        focus = "level",
    log:
        os.path.join(
            "{folder}",
            "plot_sensitivity_vs_precision",
            "logs",
            "barplot.log",
        ),
    resources:
        mem=40,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/sensitivity_vs_precision.py"

##########################################################################
##########################################################################

rule sensitivity_length:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                "{folder}",
                f"{profiler}",
                f"{size}-lognormal-{wildcards.deamination.lower()}-deamination-2.{profiler}.taxonomy.parquet", # 2 = 1x for phage
            ) for profiler in [
                "metabuli",
                "krakenuniq",
                "kraken2",
                "centrifuge",
            ] for size in ["short", "medium", "long"]
        ]
    output:
        pdf = os.path.join(
            "{folder}",
            "plot_sensitivity_newsize",
            "{deamination}_deamination.sensitivity_length.1x.pdf",
        ),
        png = os.path.join(
            "{folder}",
            "plot_sensitivity_newsize",
            "{deamination}_deamination.sensitivity_length.1x.png",
        )
    params:
        tsv = os.path.join(
            "{folder}",
            "plot_sensitivity_newsize",
            "{deamination}_deamination.sensitivity_length.parquet",
        ),
    log:
        os.path.join(
            "{folder}",
            "plot_sensitivity_newsize",
            "logs",
            "{deamination}.lineplot.log",
        ),
    resources:
        mem=80,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/sensitivity_length.py"

##########################################################################
##########################################################################

rule precision_length:
    input: 
        taxonomy = lambda wildcards: [
            os.path.join(
                "{folder}",
                f"{profiler}",
                f"{size}-lognormal-{wildcards.deamination.lower()}-deamination-2.{profiler}.taxonomy.parquet", # 2 = 1x for phage
            ) for profiler in [
                "metabuli",
                "krakenuniq",
                "kraken2",
                "centrifuge",
            ] for size in ["short", "medium", "long"]
        ]
    output:
        pdf = os.path.join(
            "{folder}",
            "plot_precision_newsize",
            "{deamination}_deamination.precision_length.1x.pdf",
        ),
        png = os.path.join(
            "{folder}",
            "plot_precision_newsize",
            "{deamination}_deamination.precision_length.1x.png",
        )
    params:
        tsv = os.path.join(
            "{folder}",
            "plot_precision_newsize",
            "{deamination}_deamination.precision_length.parquet",
        ),
    log:
        os.path.join(
            "{folder}",
            "plot_precision_newsize",
            "logs",
            "{deamination}.lineplot.log",
        ),
    resources:
        mem=80,
        cpus=2,
    threads: 2
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/precision_length.py"

##########################################################################
##########################################################################