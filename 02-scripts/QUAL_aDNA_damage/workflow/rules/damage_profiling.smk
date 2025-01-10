##########################################################################
##########################################################################

rule samtools_calmd:
    input:
        bam=os.path.join(
            BAM_FOLDER,
            "{sample}.viruses_decoy.sam2lca.total.sorted.bam",
        ),
        refs=config['refs_virus'],
    output:
        bam = os.path.join(
                OUTPUT_FOLDER,
                "BAM_files",
                "{sample}.viruses_decoy.sam2lca.total.calmd.bam",
            ),
        idx = os.path.join(
                OUTPUT_FOLDER,
                "BAM_files",
                "{sample}.viruses_decoy.sam2lca.total.calmd.bam.bai",
            ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "samtools.calmd.{sample}.log",
        ),
    resources:
        mem=10,
        cpus=5,
    threads: 5  # Use at least two threads
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        samtools calmd --threads {threads} -b {input.bam} {input.refs} > {output.bam} 2> {log:q}
        samtools index {output.bam} &>> {log:q}
        """

##########################################################################
##########################################################################

rule subset_to_mag:
    input:
        contigs_list=os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.txt",
        ),
        bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.viruses_decoy.sam2lca.total.calmd.bam",
        ),
        bai=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.viruses_decoy.sam2lca.total.calmd.bam.bai",
        ),
    output:
        bam = temp(
            os.path.join(
                OUTPUT_FOLDER,
                "BAM_files",
                "{sample}.viruses.subset.bam",
            ),
        ),
        bai = temp(
            os.path.join(
                OUTPUT_FOLDER,
                "BAM_files",
                "{sample}.viruses.subset.bam.bai",
            ),
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "samtools.subset.{sample}.log",
        ),
    resources:
        mem=4,
        cpus=5,
    threads: 5  # Use at least two threads
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        cat {input.contigs_list:q} | xargs samtools view -bh {input.bam:q} &> {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule samtools_faidx:
    input:
        os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.fasta",
        ),
    output:
        os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.fasta.fai",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.faidx.log",
        ),
    resources:
        mem=10,
        cpus=4,
    threads: 4
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule damageprofiler:
    input: 
        bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.viruses.subset.bam",
        ),
        fa=os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.fasta",
        ),
        fai=os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.fasta.fai",
        ),
        contigs_list=os.path.join(
            REFS_FOLDER,
            "{sample}",
            "{sample}.refs.txt",
        ),
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "damage_profiling",
            "{sample}",
            "DamageProfiler.log",
        )
    params:
        outdir=os.path.join(
            OUTPUT_FOLDER,
            "damage_profiling",
            "{sample}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "damage_profiling",
            "{sample}.damageprofiler.log",
        ),
    resources:
        mem=10,
        cpus=1,
    threads: 1
    conda: 
        "../envs/damageprofiler.yaml"
    shell: 
        """
        damageprofiler -i {input.bam} \
            -o {params.outdir} \
            -r {input.fa} \
            -sf {input.contigs_list} &> {log:q}
        """
##########################################################################
##########################################################################

rule pydamage:
    input: 
        bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.viruses.subset.bam",
        ),
        baifile=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.viruses.subset.bam.bai",
        ),
    output: 
        csv = os.path.join(
            OUTPUT_FOLDER,
            "pydamage",
            "{sample}",
            "pydamage_results.csv",
        ),
        filter_csv = os.path.join(
            OUTPUT_FOLDER,
            "pydamage",
            "{sample}",
            "pydamage_filtered_results.csv",
        )
    params:
        outdir=os.path.join(
            OUTPUT_FOLDER,
            "pydamage",
            "{sample}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "pydamage",
            "{sample}.pydamage.log",
        ),
    resources:
        mem=20,
        cpus=5,
    threads: 5
    conda: 
        "../envs/pydamage.yaml"
    shell: 
        """
        pydamage -o {params.outdir} analyze --force --plot {input.bam} &> {log:q}
        pydamage -o {params.outdir} filter {output.csv} &>> {log:q}
        """

##########################################################################
##########################################################################

