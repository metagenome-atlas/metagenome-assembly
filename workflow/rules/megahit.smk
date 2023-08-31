assembly_params["megahit"] = {
    "default": "",
    "meta-sensitive": "--presets meta-sensitive",
    "meta-large": " --presets meta-large",
}

if PAIRED_END and config["merge_pairs_before_assembly"]:
    ASSEMBLY_FRACTIONS = ["R1", "R2", "me"]
else:
    ASSEMBLY_FRACTIONS = MULTIFILE_FRACTIONS


localrules:
    merge_se_me_for_megahit,


rule merge_se_me_for_megahit:
    input:
        expand(
            "Intermediate/Assembly/{{sample}}/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=["se", "me"],
            assembly_preprocessing_steps=assembly_preprocessing_steps,
        ),
    output:
        temp(
            expand(
                "Intermediate/Assembly/{{sample}}/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
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
            "Intermediate/Assembly/{{sample}}/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=ASSEMBLY_FRACTIONS,
            assembly_preprocessing_steps=assembly_preprocessing_steps,
        ),
    output:
        temp("Intermediate/Assembly/megahit/{sample}/{sample}_prefilter.contigs.fa"),
    benchmark:
        "logs/benchmarks/assembly/megahit/{sample}.txt"
    log:
        "{sample}/logs/assembly/megahit.log",
    params:
        min_count=config["megahit_min_count"],
        k_min=config["megahit_k_min"],
        k_max=config["megahit_k_max"],
        k_step=config["megahit_k_step"],
        merge_level=config["megahit_merge_level"],
        prune_level=config["megahit_prune_level"],
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
        time=config["runtime_assembly"],
    shell:
        "rm -r {params.outdir} 2> {log} "
        " ;\n "
        " megahit "
        " {params.inputs} "
        " --tmp-dir {resources.tmpdir} "
        " --num-cpu-threads {threads} "
        " --k-min {params.k_min} "
        " --k-max {params.k_max} "
        " --k-step {params.k_step} "
        " --out-dir {params.outdir} "
        " --out-prefix {wildcards.sample}_prefilter "
        " --min-contig-len {params.min_contig_len} "
        " --min-count {params.min_count} "
        " --merge-level {params.merge_level} "
        " --prune-level {params.prune_level} "
        " --low-local-ratio {params.low_local_ratio} "
        " --memory {resources.mem}000000000  "
        " {params.preset} &>> {log} "


# TODO: is it necessary to remove the output dir


# standardizes header labels within contig FASTAs
rule rename_megahit_contigs:
    input:
        rules.run_megahit.output[0]
    output:
        final_assembly,
    conda:
        "../envs/bbmap.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime_simplejob"],
    log:
        "logs/assembly/post_process/rename_and_filter_size/{sample}.log",
    params:
        minlength=config["minimum_contig_length"],
    shell:
        "rename.sh "
        " in={input} out={output} ow=t "
        " prefix={wildcards.sample} "
        " minscaf={params.minlength} &> {log} "

# TODO: start with 1 not zero