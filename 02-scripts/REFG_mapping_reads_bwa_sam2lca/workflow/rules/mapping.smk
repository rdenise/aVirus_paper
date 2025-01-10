##########################################################################
##########################################################################

rule bwa_index:
    input:
        ref=os.path.join(
            REFS_FOLDER,
            "{genome}.fasta",
        )
    output:
        temp(
            multiext(
                os.path.join(
                    REFS_FOLDER,
                    "{genome}.fasta",
                ),
                ".bwt",
                ".pac",
                ".sa",
            )
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bwa",
            "bwa.index.{genome}.log",
        ),
    threads: 1  # Use at least two threads
    resources:
        mem=4,
        cpus=1,
    conda: 
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input.ref:q} &> {log:q}
        """

##########################################################################
##########################################################################

rule bwa_aln:
    input:
        sample_R1=lambda wildcards: sample2file[wildcards.sample][0],
        sample_R2=lambda wildcards: sample2file[wildcards.sample][1],
        sample_merge=lambda wildcards: sample2file[wildcards.sample][2],
        bwa_db=multiext(
            os.path.join(
                REFS_FOLDER,
                "{reference}.fasta",
            ),
                ".bwt",
                ".pac",
                ".sa",
        ),
        ref=os.path.join(
            REFS_FOLDER,
            "{reference}.fasta",
        ),
    output:
        temp(
            os.path.join(
                OUTPUT_FOLDER,
                "BAM_files",
                "{sample}.{reference}.pipe.bam",
            ),
        ),
    params:
        extra="-n 0.01 -k 2 -o 2 -l 1024",  # optional parameters
        R1_sai=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.R1.sai",
        ),
        R2_sai=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.R2.sai",
        ),
        unpaired_sai=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.unpaired.sai",
        ),
        paired_bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.paired.bam",
        ),
        unpaired_bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.unpaired.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bwa",
            "bwa.{reference}.{sample}.log",
        ),
    resources:
        mem=50,
        cpus=10,
    threads: 10  # Use at least two threads
    conda: 
        "../envs/bwa.yaml"
    shell:
        """
        # Step 1: Generate .sai alignment files for each read file
        echo "Aligning reads..." > {log:q}

        if [ -f {params.R1_sai:q} ]; then
            echo "R1 sai file exists"
        else
            bwa aln -t {threads} {params.extra} {input.ref:q} {input.sample_R1:q} > {params.R1_sai:q} 2>> {log:q}
        fi

        if [ -f {params.R2_sai:q} ]; then
            echo "R2 sai file exists"
        else
            bwa aln -t {threads} {params.extra} {input.ref:q} {input.sample_R2:q} > {params.R2_sai:q} 2>> {log:q}
        fi

        # If you have unpaired reads (single-end), align them as well
        echo "Aligning unpaired reads..." >> {log:q}

        if [ -f {params.unpaired_sai:q} ]; then
            echo "Unpaired sai file exists"
        else
            bwa aln -t {threads} {params.extra} {input.ref:q} {input.sample_merge:q} > {params.unpaired_sai:q} 2>> {log:q}
        fi

        # Step 2: Generate the SAM file
        echo "Generating SAM file..." >> {log:q}

        bwa sampe -n 50 {input.ref:q} {params.R1_sai:q} {params.R2_sai:q} {input.sample_R1:q} {input.sample_R2:q} | samtools view -Sb - > {params.paired_bam:q}  2>> {log:q}

        # For unpaired reads, use bwa samse
        echo "Generating SAM file for unpaired reads..." >> {log:q}

        bwa samse -n 50 {input.ref:q} {params.unpaired_sai:q} {input.sample_merge:q} | samtools view -Sb - > {params.unpaired_bam:q} 2>> {log:q}

        # Step 3: Merge the BAM files
        echo "Merging BAM files..." >> {log:q}

        samtools merge -@ {threads} {output:q} {params.paired_bam:q} {params.unpaired_bam:q} 2>> {log:q}

        # Step 4: Remove the intermediate files
        echo "Removing intermediate files..." >> {log:q}
        rm {params.paired_bam} {params.unpaired_bam} 
        """

##########################################################################
##########################################################################

rule samtools_view_index:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.pipe.bam",
        ),
    output:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
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
        "../envs/bwa.yaml"
    shell:
        """
        samtools view {params.extra} -b {input:q} --threads {threads} | samtools sort --threads {threads} - > {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule bwa_bam_process:
    input:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.reads.sorted.bam",
        ),
    output:
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.bam",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bwa",
            "{sample}.{reference}.bam.process.log",
        ),
    resources:
        mem=50,
        cpus=1,
    threads: 1
    wildcard_constraints:
        sample = "[^.]+",
        reference = "[^.]+",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bwa_bam_process.py"

##########################################################################
##########################################################################


rule samtools_bwa_sort_index:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.bam",
        ),
        ref = os.path.join(
            REFS_FOLDER,
            "{reference}.fasta",
        )
    output: 
        bai = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.sorted.bam.bai",
        ),
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.sorted.bam",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bamAlignCleaner",
            "{sample}.{reference}.samtools.index.log",
        ),
    wildcard_constraints:
        sample = "[^.]+",
        reference = "[^.]+",
    resources:
        mem=4,
        cpus=4,
    threads: 4
    conda: 
        "../envs/bwa.yaml"
    shell: 
        """
        samtools sort -u {input.bam:q} | samtools calmd -bAr -@ {threads} - {input.ref:q} > {output.bam:q} 2> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################


rule bamAlignCleaner:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.sorted.bam",
        ),
        bai = os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.processed.sorted.bam.bai",
        )
    output: 
        bam=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.cleaned.processed.sorted.bam",
        ),
        bai=os.path.join(
            OUTPUT_FOLDER,
            "BAM_files",
            "{sample}.{reference}.cleaned.processed.sorted.bam.bai",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bamAlignCleaner",
            "{sample}.{reference}.bamAlignCleaner.log",
        ),
    wildcard_constraints:
        sample = "[^.]+",
        reference = "[^.]+",
    resources:
        mem=4,
        cpus=4,
    threads: 4
    conda: 
        "../envs/bamAlignCleaner.yaml"
    shell: 
        """
        bamAlignCleaner {input.bam:q} | samtools sort -o {output.bam:q} &> {log:q}
        samtools index {output.bam:q} &>> {log:q}
        """

##########################################################################
##########################################################################

rule depth_breadth:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.{rank}.sorted.bam",
        ),
        fasta = config['refs_virus'],
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth.tsv",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth.log",
        ),
    resources:
        mem=4,
        cpus=1,
    threads: 1
    conda: 
        "../envs/cmseq.yaml"
    shell: 
        """
        breadth_depth.py -c {input.fasta} {input.bam} 1> {output:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule depth_breadth_cov10:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "sam2lca",
            "{sample}.{reference}.sam2lca.{rank}.sorted.bam",
        ),
        fasta = config['refs_virus'],
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth_cov10.tsv",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "depth_breadth",
            "{sample}.{reference}.{rank}.depth_breadth.cov10.log",
        ),
    resources:
        mem=4,
        cpus=1,
    threads: 1
    conda: 
        "../envs/cmseq.yaml"
    shell: 
        """
        breadth_depth.py --mincov 10 -c {input.fasta} {input.bam} 1> {output:q} 2> {log:q}
        """

##########################################################################
##########################################################################
