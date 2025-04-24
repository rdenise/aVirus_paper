##########################################################################
##########################################################################

rule uncompress:
    input: 
        refs=os.path.join(
            CONTIGS_FOLDER,
            "{genome}.fna.gz",
        )
    output: 
        fasta=temp(os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "fasta_tmp",
            "{genome}.fna"
        ))
    conda: 
        "../envs/python.yaml"
    log: 
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "uncompress",
            "{genome}.log",
        )
    resources:
        mem=10,
        cpus=1,
    threads: 1
    script:
        "../scripts/uncompress_rename.py"

##########################################################################
##########################################################################

rule concat:
    input:
        refs=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "simulation_genomad",
                "fasta_tmp",
                "{genome}.fna"
            ),
            genome=SAMPLE_NAMES,
        )
    output: 
        fasta=temp(os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "fasta_tmp",
            "allcontigs.fna",
        ))
    log: 
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "uncompress",
            "concat.log",
        )
    resources:
        mem=10,
        cpus=1,
    threads: 1
    shell:
        "cat {input.refs} 1> {output.fasta} 2> {log:q}"

##########################################################################
##########################################################################


rule genomad_default:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "fasta_tmp",
            "{genome}.fna",
        ),
        db = GENOMAD_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "genomad_default",
            "{genome}",
            "{genome}_summary",
            "{genome}_virus.fna",
        ),
    params:
        genomad_folder = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "genomad_default",
            "{genome}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "genomad_default",
            "{genome}.log",
        ),
    resources:
        mem=80,
        cpus=7,
    threads: 7
    wildcard_constraints:
        genome = "[^./]+",
        genomad = "[^./]+"
    conda: 
        "../envs/genomad.yaml"
    shell:
        """
        # genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
        genomad end-to-end --lenient-taxonomy --enable-score-calibration --composition metagenome -t {threads} --cleanup {input.contigs:q} {params.genomad_folder:q} {input.db:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule genomad:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "fasta_tmp",
            "{genome}.fna",
        ),
        db = GENOMAD_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "genomad",
            "{genome}",
            "{genome}_summary",
            "{genome}_virus.fna",
        ),
    params:
        genomad_folder = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "genomad",
            "{genome}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "genomad",
            "{genome}.log",
        ),
    resources:
        mem=80,
        cpus=7,
    threads: 7
    wildcard_constraints:
        genome = "[^./]+",
        genomad = "[^./]+"
    conda: 
        "../envs/genomad.yaml"
    shell:
        """
        # genomad end-to-end [OPTIONS] INPUT OUTPUT DATABASE
        genomad end-to-end --relaxed --lenient-taxonomy --restart --enable-score-calibration --composition metagenome -t {threads} --cleanup {input.contigs:q} {params.genomad_folder:q} {input.db:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule checkv:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "{genomad}",
            "{genome}",
            "{genome}_summary",
            "{genome}_virus.fna",
        ),
        db = CHECKV_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "checkv_{genomad}",
            "{genome}",
            "quality_summary.tsv"
        ),
    params:
        checkv_folder = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "checkv_{genomad}",
            "{genome}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "checkv_{genomad}",
            "{genome}.log",
        ),
    resources:
        mem=80,
        cpus=10,
    threads: 10
    wildcard_constraints:
        genome = "[^./]+",
        genomad = "g[^./]+"
    conda: 
        "../envs/checkv.yaml"
    shell:
        """
        # checkv end_to_end <input> <output> [options]
        checkv end_to_end {input.contigs} {params.checkv_folder} -t {threads} -d {input.db:q} --remove_tmp &> {log:q}
        """

##########################################################################
##########################################################################

rule jeager:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "fasta_tmp",
            "{genome}.fna",
        ),
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "jeager",
            "split_contigs",
            "{genome}",
            "{genome}_default_phages_jaeger.fasta",
        ),
    params:
        jeager_folder = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "jeager",
            "split_contigs",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "jeager",
            "{genome}.log",
        ),
    resources:
        mem=50,
        cpus=10,
    threads: 10
    wildcard_constraints:
        genome = "[^./]+",
    conda: 
        "../envs/jeager.yaml"
    shell:
        """
        # export OMP_NUM_THREADS={threads}
        # export TF_NUM_INTEROP_THREADS={threads}
        # export TF_NUM_INTRAOP_THREADS={threads}
        
        jaeger run -i {input.contigs} -o {params.jeager_folder} --fsize 500 --stride 500 --cpu --getalllabels --overwrite --workers {threads} --prophage --getsequences -f &> {log:q}
        """

##########################################################################
##########################################################################

rule checkv_jeager:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "jeager",
            "split_contigs",
            "{genome}",
            "{genome}_default_phages_jaeger.fasta",
        ),
        db = CHECKV_DB,
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "checkv_jeager",
            "{genome}",
            "quality_summary.tsv"
        ),
    params:
        checkv_folder = os.path.join(
            OUTPUT_FOLDER,
            "simulation_genomad",
            "checkv_jeager",
            "{genome}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "simulation_genomad",
            "checkv_jeager",
            "{genome}",
            "{genome}.log",
        ),
    resources:
        mem=80,
        cpus=10,
    threads: 10
    wildcard_constraints:
        genome = "[^./]+",
    conda: 
        "../envs/checkv.yaml"
    shell:
        """
        # checkv end_to_end <input> <output> [options]
        checkv end_to_end {input.contigs} {params.checkv_folder} -t {threads} -d {input.db:q} --remove_tmp &> {log:q}
        """

##########################################################################
##########################################################################