##########################################################################
##########################################################################

rule uncompress:
    input:
        os.path.join(
            config["metagenomes"]["assembled_contigs"],
            "{software}",
            "{sample}-{software}.fasta.gz"
        )
    output:
        temp(
            os.path.join(
                OUTPUT_FOLDER,
                "{software}",
                "{sample}-{software}.fasta"
            )
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "gzip",
            "{sample}-{software}.log",
        ),
    shell:
        "gunzip -c {input:q} > {output:q} 2> {log:q}"
            
##########################################################################
##########################################################################

rule ungzip:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna.gz"
        )
    output:
        temp(os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna"
        ))
    shell:
        "gunzip -c {input} > {output}"

##########################################################################
##########################################################################

rule vamb:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna"
        ),
        genomad=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "genomad",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna"
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "vamb",
            "{site}-{software}.fasta",
        ),
    params:
        clusters=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "vamb",
            "{site}-{software}",
            "vae_clusters_split.tsv"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "{method}",
            "vamb.concat.{software}.{site}.log",
        )
    threads: 1
    resources:
        mem=80,
        cpus=1,
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/concat_bin.py"


##########################################################################
##########################################################################

rule taxvamb_site:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna"
        ),
        genomad=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "genomad",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna"
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "taxvamb",
            "{site}-{software}.fasta",
        ),
    params:
        clusters=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "taxvamb",
            "{site}-{software}",
            "vaevae_clusters_split.tsv"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "{method}",
            "taxvamb.concat.{software}.{site}.log",
        )
    threads: 1
    wildcard_constraints:
        site="ASM|BSM"
    resources:
        mem=80,
        cpus=1,
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/concat_bin.py"

##########################################################################
##########################################################################

rule taxvamb_split_sample:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fasta"
        ),
        genomad=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "genomad",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fasta"
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "taxvamb",
            "{site}-{software}.fasta",
        ),
    params:
        clusters=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "taxvamb",
            "{site}-{software}",
            "vaevae_clusters_split.tsv"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "{method}",
            "taxvamb.concat.{software}.{site}.log",
        )
    threads: 1
    wildcard_constraints:
        site="H.*|Z.*"
    resources:
        mem=80,
        cpus=1,
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/concat_bin.py"


##########################################################################
##########################################################################

def get_path_sample_contigs(wildcards):
    if wildcards.method == "genomad":
        # return os.path.join(
        #     config["metagenomes"]["genomad_contigs"],
        #     f"{wildcards.sample}-{wildcards.software}",
        #     f"{wildcards.sample}-{wildcards.software}_summary",
        #     f"{wildcards.sample}-{wildcards.software}_virus.fna"
        # )
        return os.path.join(
                config["metagenomes"]["genomad_contigs"],
                f"{wildcards.sample}-{wildcards.software}",
                f"{wildcards.sample}-{wildcards.software}_summary",
                f"{wildcards.sample}-{wildcards.software}_virus.fna")

    elif wildcards.method == "assembled":
        # return os.path.join(
        #     config["metagenomes"]["assembled_contigs"],
        #     f"{wildcards.software}",
        #     f"{wildcards.sample}-{wildcards.software}.fasta"
        # )
        return os.path.join(
                OUTPUT_FOLDER,
                f"{wildcards.software}",
                f"{wildcards.sample}-{wildcards.software}.fasta"
            )

##########################################################################

def get_genomad_fasta(wildcards):
    return os.path.join(
        config["metagenomes"]["genomad_contigs"],
        f"{wildcards.sample}-{wildcards.software}",
        f"{wildcards.sample}-{wildcards.software}_summary",
        f"{wildcards.sample}-{wildcards.software}_virus.fna"
    )

##########################################################################
##########################################################################

rule semibin:
    input:
        contigs = lambda wildcards: get_path_sample_contigs(wildcards),
        genomad = lambda wildcards: get_genomad_fasta(wildcards),
    output:
        fasta = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "semibin",
            "{sample}-{software}.fasta",
        )
    params:
        clusters = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "semibin",
            "{sample}-{software}",
            "contig_bins.tsv"
        ) 
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "semibin",
            "{sample}-{software}.log",
        )
    threads: 1
    resources:
        mem=80,
        cpus=1,
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/concat_bin.py"

##########################################################################
##########################################################################

rule split:
    input:
        contigs = lambda wildcards: os.path.join(
            OUTPUT_FOLDER,
            "binning",
            f"{wildcards.method}",
            f"{wildcards.software}",
            "concat_bins_fasta",
            f"{wildcards.profiler}",
            f"{wildcards.sample[:3]}-{wildcards.software}.fasta",
        )
    output:
        fasta = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "split",
            "{profiler}",
            "{sample}-{software}.log",
        )
    threads: 1
    resources:
        mem=10,
        cpus=1,
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/split_vamb.py"

##########################################################################
##########################################################################

rule checkv:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        ),
        db = config["checkv_db"],
    output:
        summary = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "checkv",
            "{sample}-{software}",
            "quality_summary.tsv"
        ),
    params:
        checkv_folder = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "checkv",
            "{sample}-{software}",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "{profiler}",
            "checkv",
            "{sample}.checkv.{method}.{software}.log",
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
        if [ -s {input.contigs} ]; then
            checkv end_to_end {input.contigs} {params.checkv_folder} -t {threads} -d {input.db:q} --remove_tmp &> {log:q}
        else
            mkdir -p {params.checkv_folder}
            touch {output.summary}
        fi
        """

##########################################################################
##########################################################################

rule bwa_index:
    input:
        ref=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        ),
    output:
        temp(
            multiext(
                os.path.join(
                    OUTPUT_FOLDER,
                    "binning",
                    "{method}",
                    "{software}",
                    "concat_bins_fasta",
                    "{profiler}",
                    "{sample}-{software}.fasta",
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
            "{profiler}",
            "bwa.index.{method}.{sample}-{software}.log",
        ),
    threads: 1 
    resources:
        mem=4,
        cpus=1,
    conda: 
        "../envs/bwa.yaml"
    shell:
        """
        if [ -s {input.ref} ]; then
            bwa index {input.ref:q} &> {log:q}
        else
            touch {input.ref:q}.bwt {input.ref:q}.pac {input.ref:q}.sa
        fi
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
                "binning",
                "{method}",
                "{software}",
                "concat_bins_fasta",
                "{profiler}",
                "{sample}-{software}.fasta",
            ),
            ".bwt",
            ".pac",
            ".sa",
        ),
        ref=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        ),
    output:
        temp(
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                "{method}",
                "{software}",
                "concat_bins_fasta",
                "{profiler}",
                "pydamage",
                "BAM_files",
                "{sample}-{software}.pipe.bam",
            ),
        ),
    params:
        extra="-n 0.01 -k 2 -o 2 -l 1024",  # optional parameters
        R1_sai=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.R1.sai",
        ),
        R2_sai=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.R2.sai",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "{profiler}",
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
        if [ -s {input.ref} ]; then
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
            bwa sampe {input.ref:q} {params.R1_sai:q} {params.R2_sai:q} {input.sample_R1:q} {input.sample_R2:q} | samtools view -Sb - > {output:q} 2>> {log:q}
        else
            touch {output:q}
        fi
        """

##########################################################################
##########################################################################


rule samtools_bwa_sort_index:
    input: 
        bam = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.pipe.bam",
        ),
        ref = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        ),
    output: 
        bai = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam.bai",
        ),
        bam = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
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
            "{profiler}",
            "{method}",
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
        if [ -s {input.ref} ]; then
            samtools sort -u {input.bam:q} | samtools calmd -bAr -@ {threads} - {input.ref:q} > {output.bam:q} 2> {log:q}
            samtools index {output.bam:q} &>> {log:q}
        else
            touch {output.bam:q} {output.bai:q}
        fi
        """

##########################################################################
##########################################################################

rule pydamage:
    input: 
        bam=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam",
        ),
        baifile=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "BAM_files",
            "{sample}-{software}.calmd.sorted.bam.bai",
        ),
    output: 
        csv = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "{sample}-{software}",
            "pydamage_results.csv",
        ),
        filter_csv = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "{sample}-{software}",
            "pydamage_filtered_results.csv",
        )
    params:
        outdir=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "pydamage",
            "{sample}-{software}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "{profiler}",
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
        if [ -s {input.bam} ]; then
            pydamage -o {params.outdir} analyze -p {threads} --force {input.bam} &> {log:q}
            pydamage -o {params.outdir} filter {output.csv} &>> {log:q}
        else
            mkdir -p {params.outdir}
            touch {output.csv} {output.filter_csv}
        fi
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
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "{sample}-{software}.fasta",
        )
    output:
        ani = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.ani.tsv",
        ),
        ani_ids = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.ani.ids.tsv",
        ),
        filter_out = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.filter.tsv",
        ),
        clusters = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.clusters.tsv",
        ),
        aln = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.contigs.{software}.aln.tsv",
        )
    params:
        tmp = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concat_bins_fasta",
            "{profiler}",
            "vclust",
            "{sample}.{software}.tmp",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "{software}",
            "vclust",
            "{sample}.{profiler}.vclust.{software}.log",
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
        if [ -s {input.query} ]; then
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
        else
            touch {output:q}
        fi
        """

##########################################################################
##########################################################################