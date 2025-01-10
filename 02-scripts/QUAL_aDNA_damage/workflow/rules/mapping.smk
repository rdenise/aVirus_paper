##########################################################################
##########################################################################

rule bowtie2_index:
    input:
        ref=os.path.join(
            REFS_FOLDER_VIRUS,
            "{genome}.fasta",
        )
    output:
        temp(
            multiext(
                os.path.join(
                    REFS_FOLDER_VIRUS,
                    "{genome}.fasta",
                ),
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            )
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2",
            "bowtie2.index.{genome}.log",
        ),
    threads: 5  # Use at least two threads
    conda: 
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --threads {threads} {input.ref:q} {input.ref:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule bowtie2:
    input:
        sample_R1=lambda wildcards: sample2file[wildcards.sample][0],
        sample_R2=lambda wildcards: sample2file[wildcards.sample][1],
        bwt_db=multiext(
            os.path.join(
                REFS_FOLDER_VIRUS,
                "{reference}.fasta",
            ),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        ref=os.path.join(
            REFS_FOLDER_VIRUS,
            "{reference}.fasta",
        ),
    output:
        temp(
            os.path.join(
                BAM_FOLDER,
                "{sample}.{reference}.pipe.bam",
            ),
        ),
    params:
        # extra="-D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --no-unal",  # optional parameters
        extra="-D 20 -R 3 -N 1 -L 20 -k 50 -i S,1,0.50 --no-unal",  # optional parameters
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2",
            "bowtie2.{reference}.{sample}.log",
        ),
    resources:
        mem=4,
        cpus=10,
    threads: 10  # Use at least two threads
    conda: 
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {input.ref:q} -1 {input.sample_R1:q} -2 {input.sample_R2:q} {params.extra} -p {threads} | samtools view -Sb - > {output:q} 2> {log:q}"

##########################################################################
##########################################################################

rule samtools_view_index:
    input:
        os.path.join(
            BAM_FOLDER,
            "{sample}.{reference}.pipe.bam",
        ),
    output:
        bam = os.path.join(
            BAM_FOLDER,
            "{sample}.{reference}.reads.sorted.bam",
        ),
    params:
        extra="-F 4",  # optional params string
        region="",  # optional region string
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
        "../envs/bowtie2.yaml"
    shell:
        """
        samtools view {params.extra} -b {input:q} --threads {threads} | samtools sort --threads {threads} - > {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule bamAlignCleaner:
    input: 
        os.path.join(
            BAM_FOLDER,
            "{sample}.{reference}.reads.sorted.bam",
        ),
    output: 
        os.path.join(
            BAM_FOLDER,
            "{sample}.{reference}.cleaned.sorted.bam",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bamAlignCleaner",
            "{sample}.{reference}.bamAlignCleaner.log",
        ),
    resources:
        mem=4,
        cpus=4,
    threads: 4
    conda: 
        "../envs/bamAlignCleaner.yaml"
    shell: 
        """
        bamAlignCleaner {input:q} | samtools sort -o {output:q} &> {log:q}
        samtools index {output:q} &>> {log:q}
        """

##########################################################################
##########################################################################