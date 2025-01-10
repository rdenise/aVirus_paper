##########################################################################
##########################################################################

rule sam2lca:
    input:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.cleaned.processed.sorted.bam",
        ),
        bai = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.cleaned.processed.sorted.bam.bai",
        ),
    output:
        csv = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.csv",
        ),
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.bam",
        ),
    params:
        sam2lca_db = SAM2LCA_DB,
        seqid2taxid = SEQID2TAXID,
        taxonomy = TAXONOMY,
        sam2lca_folder = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca"
        ),
        sample = '{sample}',
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "sam2lca",
            "{sample}.{reference}.sam2lca.analyze.log",
        ),
    resources:
        mem=120,
        cpus=4,
    threads: 4  
    conda:
        "../envs/sam2lca.yaml"
    shell:
        """
        sam2lca -d {params.sam2lca_db:q} analyze \
        -p {threads} -b \
        -i 0.9 \
        -o {output.csv:q} \
        --acc2tax {params.seqid2taxid} \
        --taxonomy {params.taxonomy} \
        {input.bam:q} &> {log:q}

        rename 's/sam2lca.sam2lca/sam2lca/' * &>> {log:q}
        mv {params.sample}*sam2lca* {params.sam2lca_folder:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule split_sam2lca_BAM:
    input:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.bam",
        )
    output:
        genus = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.genus.bam",
        ),
        other = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.higher.bam",
        ),
        species = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.species.bam",
        ),
        unclassified = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.unclassified.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "sam2lca",
            "{sample}.{reference}.bam.split.log",
        ),
    resources:
        mem=20,
        cpus=10,
    threads: 10  
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/postprocess_bam.py"

##########################################################################
##########################################################################

rule samtools_sort_index:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.{rank}.bam",
        ),
    output:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.{rank}.sorted.bam",
        ),
    params:
        extra="-F 4",  # optional params string
        region="",  # optional region string
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.{reference}.{rank}.samtools.view.log",
        ),
    resources:
        mem=10,
        cpus=4,
    threads: 4
    conda: 
        "../envs/bwa.yaml"
    shell:
        """
        samtools view {params.extra} -b {input:q} --threads {threads} | samtools sort --threads {threads} - > {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule samtools_sort_index_total:
    input:
        genus = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.genus.sorted.bam",
        ),
        other = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.higher.sorted.bam",
        ),
        species = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.species.sorted.bam",
        ),
    output:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.total.sorted.bam",
        ),
    params:
        extra="-F 4",  # optional params string
        region="",  # optional region string
        tmp_merge = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.tmp.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.{reference}.samtools.view.log",
        ),
    resources:
        mem=10,
        cpus=4,
    threads: 4
    conda: 
        "../envs/bwa.yaml"
    shell:
        """
        echo "Merging BAM files" &> {log:q}
        samtools merge -f {params.tmp_merge:q} {input.genus:q} {input.other:q} {input.species:q} &>> {log:q}

        echo "Sorting and indexing BAM file" &>> {log:q}
        samtools view {params.extra} -b {params.tmp_merge:q} --threads {threads} | samtools sort --threads {threads} - > {output.bam:q} 2>> {log:q}
        
        echo "Indexing BAM file" &>> {log:q}
        samtools index {output.bam:q} &>> {log:q}

        echo "Removing temporary files" &>> {log:q}
        rm {params.tmp_merge:q} &>> {log:q}
        """

##########################################################################
##########################################################################