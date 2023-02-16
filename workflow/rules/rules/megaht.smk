    assembly_params["megahit"] = {
        "default": "",
        "meta-sensitive": "--presets meta-sensitive",
        "meta-large": " --presets meta-large",
    }
    ASSEMBLY_FRACTIONS = MULTIFILE_FRACTIONS
    if PAIRED_END and config.get("merge_pairs_before_assembly", True):
        ASSEMBLY_FRACTIONS = ["R1", "R2", "me"]

    localrules:
        merge_se_me_for_megahit,

    rule merge_se_me_for_megahit:
        input:
            expand(
                "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=["se", "me"],
                assembly_preprocessing_steps=assembly_preprocessing_steps,
            ),
        output:
            temp(
                expand(
                    "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                    fraction=["co"],
                    assembly_preprocessing_steps=assembly_preprocessing_steps,
                )
            ),
        shell:
            "cat {input} > {output}"


    def megahit_input_parsing(input):
        Nfiles = len(input)

        if Nfiles == 1:
            out = f"--read {input[0]}"
        else:
            out = f"-1 {input[0]} -2 {input[1]} "

            if Nfiles == 3:
                out += f"--read {input[2]}"
        return out

    rule run_megahit:
        input:
            expand(
                "{{sample}}/assembly/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
                fraction=ASSEMBLY_FRACTIONS,
                assembly_preprocessing_steps=assembly_preprocessing_steps,
            ),
        output:
            temp("{sample}/assembly/megahit/{sample}_prefilter.contigs.fa"),
        benchmark:
            "logs/benchmarks/assembly/megahit/{sample}.txt"
        log:
            "{sample}/logs/assembly/megahit.log",
        params:
            min_count=config.get("megahit_min_count", MEGAHIT_MIN_COUNT),
            k_min=config.get("megahit_k_min", MEGAHIT_K_MIN),
            k_max=config.get("megahit_k_max", MEGAHIT_K_MAX),
            k_step=config.get("megahit_k_step", MEGAHIT_K_STEP),
            merge_level=config.get("megahit_merge_level", MEGAHIT_MERGE_LEVEL),
            prune_level=config.get("megahit_prune_level", MEGAHIT_PRUNE_LEVEL),
            low_local_ratio=config["megahit_low_local_ratio"],
            min_contig_len=config["minimum_contig_length"],
            outdir=lambda wc, output: os.path.dirname(output[0]),
            inputs=lambda wc, input: megahit_input_parsing(input),
            preset=assembly_params["megahit"][config["megahit_preset"]],
        conda:
            "../envs/megahit.yaml"
        threads: config["assembly_threads"]
        resources:
            mem=config["assembly_memory"],
            time=config["runtime"]["assembly"],
        shell:
            """
            rm -r {params.outdir} 2> {log}

            megahit \
            {params.inputs} \
            --tmp-dir {TMPDIR} \
            --num-cpu-threads {threads} \
            --k-min {params.k_min} \
            --k-max {params.k_max} \
            --k-step {params.k_step} \
            --out-dir {params.outdir} \
            --out-prefix {wildcards.sample}_prefilter \
            --min-contig-len {params.min_contig_len} \
            --min-count {params.min_count} \
            --merge-level {params.merge_level} \
            --prune-level {params.prune_level} \
            --low-local-ratio {params.low_local_ratio} \
            --memory {resources.mem}000000000  \
            {params.preset} >> {log} 2>&1
            """


