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
        mem=120,
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
        genomad end-to-end --relaxed --restart --lenient-taxonomy --splits 8 --enable-score-calibration --composition metagenome -t {threads} --cleanup {input.contigs:q} {params.genomad_folder:q} {input.db:q} &> {log:q}
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
            "checkv_genomad",
            "{sample}-{software}",
            "quality_summary.tsv"
        ),
    params:
        checkv_folder = os.path.join(
            OUTPUT_FOLDER,
            "checkv_genomad",
            "{sample}-{software}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv_genomad",
            "{sample}.checkv.{software}.log",
        ),
    resources:
        mem=40,
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

rule vclust:
    input:
        target=os.path.join(
            REFS_FILE
        ), 
        query=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        )
    output:
        ani = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        ani_ids = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "vclust",
            "{sample}.contigs.{software}.ani.ids.tsv",
        ),
        filter_out = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "vclust",
            "{sample}.contigs.{software}.filter.tsv",
        ),
        clusters = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "vclust",
            "{sample}.contigs.{software}.clusters.tsv",
        ),
        aln = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "vclust",
            "{sample}.contigs.{software}.aln.tsv",
        )
    params:
        tmp = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "vclust",
            "{sample}.{software}.tmp",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "genomad",
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
        cat {input.query:q} {input.target:q} > {params.tmp:q} 2>> {log:q}

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
            vclust cluster -i {output.ani:q} -o {output.clusters:q} --ids {output.ani_ids:q} --algorithm complete --metric gani --gani 0.90 &>> {log:q}
        fi

        # Remove the temporary file.
        rm {params.tmp}* &>> {log:q}
        """

##########################################################################
##########################################################################

rule bwa_index:
    input:
        ref=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        ),
    output:
        temp(
            multiext(
                os.path.join(
                OUTPUT_FOLDER,
                "genomad",
                "{sample}-{software}",
                "{sample}-{software}_summary",
                "{sample}-{software}_virus.fna",
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
            "bwa.index.{sample}-{software}.log",
        ),
    threads: 1 
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
        sample_R1=os.path.join(
            FASTQ_FOLDER,
            "{sample}_1.fastq.gz",
        ),
        sample_R2=os.path.join(
            FASTQ_FOLDER,
            "{sample}_2.fastq.gz",
        ),
        bwa_db=multiext(
            os.path.join(
                OUTPUT_FOLDER,
                "genomad",
                "{sample}-{software}",
                "{sample}-{software}_summary",
                "{sample}-{software}_virus.fna",
            ),
                ".bwt",
                ".pac",
                ".sa",
        ),
        ref=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        ),
    output:
        temp(
            os.path.join(
                OUTPUT_FOLDER,
                "genomad",
                "pydamage",
                "BAM_files",
                "{sample}-{software}.pipe.bam",
            ),
        ),
    params:
        extra="-n 0.01 -k 2 -o 2 -l 1024",  # optional parameters
        R1_sai=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.R1.sai",
        ),
        R2_sai=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.R2.sai",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "pydamage",
            "bwa",
            "bwa.{sample}-{software}.log",
        ),
    resources:
        mem=50,
        cpus=20,
    threads: 20  # Use at least two threads
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

        # Step 2: Generate the SAM file
        echo "Generating SAM file..." >> {log:q}

        bwa sampe {input.ref:q} {params.R1_sai:q} {params.R2_sai:q} {input.sample_R1:q} {input.sample_R2:q} | samtools view -Sb - > {output:q}  2>> {log:q}
        """

##########################################################################
##########################################################################


rule samtools_bwa_sort_index:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.pipe.bam",
        ),
        ref = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "{sample}-{software}",
            "{sample}-{software}_summary",
            "{sample}-{software}_virus.fna",
        ),
    output: 
        bai = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam.bai",
        ),
        bam = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "pydamge",
            "camlmd",
            "{sample}-{software}.samtools.index.log",
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

rule pydamage:
    input: 
        bam=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam",
        ),
        baifile=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam.bai",
        ),
    output: 
        csv = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "{sample}-{software}",
            "pydamage_results.csv",
        ),
        filter_csv = os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "{sample}-{software}",
            "pydamage_filtered_results.csv",
        )
    params:
        outdir=os.path.join(
            OUTPUT_FOLDER,
            "genomad",
            "pydamage",
            "{sample}-{software}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "pydamage",
            "{sample}-{software}.pydamage.log",
        ),
    resources:
        mem=20,
        cpus=20,
    threads: 20
    conda: 
        "../envs/pydamage.yaml"
    shell: 
        """
        pydamage -o {params.outdir} analyze -p {threads} --force {input.bam} &> {log:q}
        pydamage -o {params.outdir} filter {output.csv} &>> {log:q}
        """

##########################################################################
##########################################################################