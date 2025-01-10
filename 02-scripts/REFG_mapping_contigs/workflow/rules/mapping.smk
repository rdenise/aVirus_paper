##########################################################################
##########################################################################

rule vclust:
    input:
        target=os.path.join(
            REFS_FILE
        ), 
        query=os.path.join(
            CONTIGS_FOLDER,
            "{software}",
            "{sample}-{software}.fasta.gz",
        )
    output:
        ani = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        ani_ids = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.ani.ids.tsv",
        ),
        filter_out = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.filter.tsv",
        ),
        clusters = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.clusters.tsv",
        ),
        aln = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.aln.tsv",
        )
    params:
        tmp = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.{software}.tmp",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vclust",
            "{sample}.vclust.{software}.log",
        ),
    resources:
        mem=80,
        cpus=10,
    threads: 10
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    conda: 
        "../envs/vclust.yaml"
    shell:
        """
        # Concate the reference and query genomes.
        echo "Concating the reference and query genomes" &> {log:q}
        zcat {input.query:q} > {params.tmp:q} 2>> {log:q}
        cat {input.target:q} >> {params.tmp:q} 2>> {log:q}

        # Create a pre-alignment filter.
        echo "Creating a pre-alignment filter" &>> {log:q}
        if [ ! -f {output.filter_out:q} ]; then
            vclust prefilter -t {threads} -i {params.tmp:q} -o {output.filter_out:q} --min-ident 0.70 &>> {log:q}
        fi

        # Output only genome pairs with ANI ≥ 0.90 and query coverage ≥ 0.85.
        echo "Output only genome pairs with ANI ≥ 0.90 and query coverage ≥ 0.85" &>> {log:q}

        if [ ! -f {output.ani:q} ]; then
            vclust align -t {threads} -i {params.tmp:q} -o {output.ani:q} --filter {output.filter_out:q} --out-ani 0.90 --out-qcov 0.70 --outfmt complete --out-aln {output.aln} &>> {log:q}
            # vclust align -t {threads} -i {params.tmp:q} -o {output.ani:q} --out-ani 0.90 --out-qcov 0.85 --outfmt complete --out-aln {output.aln} &>> {log:q}
        fi
        
        # Cluster contigs into vOTUs using the MIUVIG standards and the Leiden algorithm.
        echo "Clustering contigs into vOTUs using the MIUVIG standards and the Leiden algorithm" &>> {log:q}
        if [ ! -f {output.clusters:q} ]; then
            vclust cluster -i {output.ani:q} -o {output.clusters:q} --ids {output.ani_ids:q} --algorithm leiden --metric ani --ani 0.90 --qcov 0.85 --leiden-resolution 0.1 &>> {log:q}
        fi

        # Remove the temporary file.
        rm {params.tmp}* &>> {log:q}
        """

##########################################################################
##########################################################################

rule add_align:
    input: 
        aln = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.aln.tsv",
        ),
        metadata = METADATA,
        contigs = os.path.join(
            CONTIGS_FOLDER,
            "{software}",
            "{sample}-{software}.fasta.gz",
        ),
        refs = os.path.join(
            REFS_FILE
        )
    output:
        parquet = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.viruses.parquet",
        ),
    resources:
        mem=20,
        cpus=10,
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vclust",
            "{sample}.vclust.{software}.add_align.log",
        ),
    conda: 
        "../envs/python.yaml"
    threads: 10
    script:
        "../scripts/add_align.py"

##########################################################################
##########################################################################

rule vclust2bam:
    input: 
        aln = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.viruses.parquet",
        ),
        ani = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        metadata = METADATA,
        fasta = os.path.join(
            CONTIGS_FOLDER,
            "{software}",
            "{sample}-{software}.fasta.gz",
        )
    output:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.{metrics}.sorted.bam",
        ),
    params:
        min_ani = 0.90,
        metrics = "{metrics}",
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    resources:
        mem=20,
        cpus=1,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "vclust",
            "{sample}.vclust.{software}.{metrics}.bam.log",
        ),
    conda: 
        "../envs/python.yaml"
    threads: 1
    script:
        "../scripts/vclust2bam.py"

##########################################################################
##########################################################################

rule samtools_calmd_sort_index:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.sorted.bam",
        ),
        ref = os.path.join(
            REFS_FILE
        ),
    output:
        bai = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.calmd.sorted.bam.bai",
        ),
        bam = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.calmd.sorted.bam",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.{software}.samtools.calmd.log",
        ),
    wildcard_constraints:
        sample = "[^./]+",
        software = "[^./]+"
    resources:
        mem=20,
        cpus=4,
    threads: 4
    conda: 
        "../envs/samtools.yaml"
    shell: 
        """
        samtools calmd -bAr -@ {threads} {input.bam:q} {input.ref:q} > {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule depth_breadth:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "{sample}.contigs.{software}.{metrics}.sorted.bam",
        ),
        fasta = os.path.join(
            REFS_FILE
        ),
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "vclust",
            "depth_breadth",
            "{sample}.{software}.{rank}.{metrics}.depth_breadth.tsv",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "depth_breadth",
            "{sample}.{software}.{rank}.{metrics}.depth_breadth.log",
        ),
    resources:
        mem=4,
        cpus=1,
    threads: 1
    conda: 
        "../envs/cmseq.yaml"
    shell: 
        """
        if samtools view -c {input.bam:q} | grep -q '^0$'; then
            echo "Contig\tBreadth\tDepth avg\tDepth median\n" 1> {output:q} 2> {log:q}
        else
            breadth_depth.py {input.bam} 1> {output:q} 2> {log:q}
        fi
        """

##########################################################################
##########################################################################