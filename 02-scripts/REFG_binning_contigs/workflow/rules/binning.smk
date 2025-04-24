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

def get_path_contigs(wildcards):
    if wildcards.method == "genomad":
        # return os.path.join(
        #     config["metagenomes"]["genomad_contigs"],
        #     f"{wildcards.sample}-{wildcards.software}",
        #     f"{wildcards.sample}-{wildcards.software}_summary",
        #     f"{wildcards.sample}-{wildcards.software}_virus.fna"
        # )
        return sorted([
            os.path.join(
                config["metagenomes"]["genomad_contigs"],
                f"{sample}-{wildcards.software}",
                f"{sample}-{wildcards.software}_summary",
                f"{sample}-{wildcards.software}_virus.fna"
            ) for sample in SAMPLES
        ])
    elif wildcards.method == "assembled":
        # return os.path.join(
        #     config["metagenomes"]["assembled_contigs"],
        #     f"{wildcards.software}",
        #     f"{wildcards.sample}-{wildcards.software}.fasta"
        # )
        return sorted([
            os.path.join(
                OUTPUT_FOLDER,
                f"{wildcards.software}",
                f"{sample}-{wildcards.software}.fasta"
            ) for sample in SAMPLES
        ])

##########################################################################
##########################################################################

rule concatenation:
    input:
        refs=lambda wildcards: get_path_contigs(wildcards),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concatenation",
            "{software}.concatenate.fasta.gz"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "concatenation",
            "{software}.concatenate.log",
        )
    threads: 1
    resources:
        mem=10,
        cpus=1,
    conda:
        "../envs/vamb.yaml"
    script: 
        "../scripts/concatenate.py"


##########################################################################
##########################################################################

rule strobealign_vamb:
    input: 
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concatenation",
            "{software}.concatenate.fasta.gz"
        ),
        sample_R1=os.path.join(
            FASTQ_FOLDER,
            "{sample}_1.fastq.gz",
        ),
        sample_R2=os.path.join(
            FASTQ_FOLDER,
            "{sample}_2.fastq.gz",
        ),
    output: 
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{sample}-{software}.abundance.tsv",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "strobealign",
            "vamb",
            "{sample}-{software}.log",
        )
    threads: 40
    resources:
        mem=120,
        cpus=40,
    conda:
        "../envs/strobealign_semibin2.yaml"
    shell: 
        "strobealign -t {threads} -R 6 --aemb {input.contigs} {input.sample_R1} {input.sample_R2} > {output.abundance} 2> {log}"

##########################################################################
##########################################################################

rule concatenate_abundances_vamb:
    input: 
        abundances=sorted(expand(
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                "{{method}}",
                "{{software}}",
                "strobealign",
                "vamb",
                "{sample}-{{software}}.abundance.tsv",
            ),
            sample=SAMPLES,
        )), 
    output: 
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{software}.abundance.tsv"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "strobealign",
            "abundance",
            "{software}.log",
        )
    resources:
        mem=10,
        cpus=1,
    threads: 1
    shell: 
        """
        HEADER="contigname"
        FILES=({input.abundances})

        for i in $(seq 1 $((${{#FILES[@]}} - 1))); do
            HEADER+="\t$i"
        done

        CMD="paste ${{FILES[0]}}"
        for file in "${{FILES[@]:1}}"; do
            CMD+=" <(cut -f 2 \"$file\")"
        done

        echo -e "$HEADER" > {output.abundance:q}
        eval "$CMD" >> {output.abundance:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule split_abundance_vamb:
    input:
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{software}.abundance.tsv"
        ),
        concatenated=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "concatenation",
            "{software}.concatenate.fasta.gz"
        )
    output:
        abundances=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                "{{method}}",
                "{{software}}",
                "strobealign",
                "vamb",
                "{sites}-{{software}}.abundance.tsv"
            ),
            sites=['ASM', 'BSM', 'HSM', 'ZSM'],
        ),
        contigs=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                "{{method}}",
                "{{software}}",
                "strobealign",
                "vamb",
                "{sites}-{{software}}.fna.gz"
            ),
            sites=['ASM', 'BSM', 'HSM', 'ZSM'],
        )
    params:
        samples=lambda wildcards: get_path_contigs(wildcards)
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "strobealign",
            "split",
            "{software}.log",
        )
    resources:
        mem=10,
        cpus=1,
    threads: 1
    conda:
        "../envs/python.yaml"
    script: 
        "../scripts/split_abundance.py"

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
            "{site}-{software}.fna.gz"
        ),
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.abundance.tsv"
        )
    output:
        clusters=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "vamb",
            "{site}-{software}",
            "vae_clusters_metadata.tsv"
        ),
    params:
        output=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "vamb",
            "{site}-{software}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "{method}",
            "vamb.{software}.{site}.log",
        )
    threads: 20
    resources:
        mem=80,
        cpus=20,
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {params.output:q}

        vamb bin default --outdir {params.output:q}\
        --fasta {input.contigs:q}\
        -m 400 -p {threads}\
        --abundance_tsv {input.abundance:q} &> {log:q}
        """


##########################################################################
##########################################################################

rule taxvamb:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.fna.gz"
        ),
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "vamb",
            "{site}-{software}.abundance.tsv"
        ),
        taxonomy=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "mmseqs",
            "{site}.mmseqs2.tsv"
        ),
    output:
        clusters=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "taxvamb",
            "{site}-{software}",
            "vaevae_clusters_metadata.tsv"
        ),
    params:
        output=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "taxvamb",
            "{site}-{software}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "{method}",
            "taxvamb.{software}.{site}.log",
        )
    threads: 20
    resources:
        mem=80,
        cpus=20,
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {params.output:q}

        vamb bin taxvamb --outdir {params.output:q}\
        --fasta {input.contigs:q}\
        --taxonomy {input.taxonomy:q}\
        -m 400 -p {threads}\
        --abundance_tsv {input.abundance:q} &> {log:q}
        """



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
##########################################################################
rule strobealign_semibin2:
    input: 
        contigs=lambda wildcards: get_path_sample_contigs(wildcards),
        sample_R1=os.path.join(
            FASTQ_FOLDER,
            "{sample2}_1.fastq.gz",
        ),
        sample_R2=os.path.join(
            FASTQ_FOLDER,
            "{sample2}_2.fastq.gz",
        ),
    output: 
        abundance=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "semibin",
            "{sample}-{software}",
            "{sample}-{software}-{sample2}.abundance.tsv",
        ),
    params:
        split_contigs=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "semibin",
            "{sample}-{software}",
            "split_contigs.fna.gz"
        ),
        outdir = os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "strobealign",
            "semibin",
            "{sample}-{software}",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "binning",
            "{method}",
            "strobealign",
            "semibin",
            "{sample}-{software}.{sample2}.log",
        )
    threads: 40
    resources:
        mem=120,
        cpus=40,
    conda:
        "../envs/strobealign_semibin2.yaml"
    shell: 
        """
        mkdir -p {params.outdir:q}

        SemiBin2 split_contigs -i {input.contigs} -o {params.outdir:q} &> {log:q}

        strobealign -t {threads} -R 6 --aemb {params.split_contigs} {input.sample_R1} {input.sample_R2} > {output.abundance} 2> {log}
        """

##########################################################################
##########################################################################

rule semibin:
    input:
        abundances=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "binning",
                "{{method}}",
                "{{software}}",
                "strobealign",
                "semibin",
                "{{sample}}-{{software}}",
                "{{sample}}-{{software}}-{sample2}.abundance.tsv",
            ),
            sample2=SAMPLES,
        ),
        contigs=lambda wildcards: get_path_sample_contigs(wildcards),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "semibin",
            "{sample}-{software}",
            "contig_bins.tsv"
        )        
    params:
        outdir=os.path.join(
            OUTPUT_FOLDER,
            "binning",
            "{method}",
            "{software}",
            "semibin",
            "{sample}-{software}",
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
    threads: 20
    conda:
        "../envs/strobealign_semibin2.yaml"
    shell:
        """
        SemiBin2 single_easy_bin --input-fasta {input.contigs:q} \
        --abundance {input.abundances:q} \
        --output {params.outdir:q} \
        -t {threads} &> {log:q}
        """

##########################################################################
##########################################################################
