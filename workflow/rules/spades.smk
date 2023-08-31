SPADES_DIR= "Intermediate/Assembly/spades/{sample}"


from copy import deepcopy

if PAIRED_END:
    ASSEMBLY_FRACTIONS = ["R1", "R2"]
    if config["merge_pairs_before_assembly"]:
        ASSEMBLY_FRACTIONS += ["me"]
else:
    ASSEMBLY_FRACTIONS = deepcopy(MULTIFILE_FRACTIONS)

    if config["spades_preset"] == "meta":
        logger.error(
            "Metaspades cannot handle single end libraries. Use another assembler or specify 'spades_preset': normal"
        )
        exit(1)

assembly_params["spades"] = {"meta": "--meta", "normal": "", "rna": "--rna"}


def spades_input(wc,input):
    if (Path(SPADES_DIR.format(**wc))/"params.txt").exists():

        return " --restart-from last "
    
    else:

        


def spades_parameters(wc, input):
    if (Path(SPADES_DIR.format(**wc))/"params.txt").exists():

        params = {
            "inputs": "--restart-from last",
            "preset": "",
            "skip_error_correction": "",
            "extra": "",
            "longreads": "",
        }

    else:
        params = {}

        # intputs
        reads = dict(zip(ASSEMBLY_FRACTIONS, input))

        if not PAIRED_END:
            params["inputs"] = " -s {se} ".format(**reads)
        else:
            params["inputs"] = " --pe1-1 {R1} --pe1-2 {R2} ".format(**reads)

            if "se" in ASSEMBLY_FRACTIONS:
                params["inputs"] += "--pe1-s {se} ".format(**reads)
            if "me" in ASSEMBLY_FRACTIONS:
                params["inputs"] += "--pe1-m {me} ".format(**reads)

        # Long reads:

        if (config["longread_type"] is not None) & (
            str(config["longread_type"]).lower() != "none"
        ):
            long_read_file = get_long_reads(wc)
            params["longreads"] = " --{t} {f} ".format(
                t=config["longread_type"], f=long_read_file
            )
        else:
            params["longreads"] = ""

        params["preset"] = assembly_params["spades"][config["spades_preset"]]
        params["skip_error_correction"] = (
            "--only-assembler" if config["spades_skip_BayesHammer"] else ""
        )
        params["extra"] = config["spades_extra"]


    return params


rule run_spades:
    input:
        expand(
            "Intermediate/Assembly/{{sample}}/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=ASSEMBLY_FRACTIONS,
            assembly_preprocessing_steps=assembly_preprocessing_steps,
        ),
    output:

        multiext("Intermediate/Assembly/spades/{sample}/", "contigs.fasta", "scaffolds.fasta",
        "assembly_graph_with_scaffolds.gfa",
        "scaffolds.paths",
        "contig.paths")
    benchmark:
        "logs/benchmarks/assembly/spades/{sample}.txt"
    params:
        p=lambda wc, input: spades_parameters(wc, input),
        outdir = SPADES_DIR,
        k=config["spades_k"],
    log:
        "{sample}/logs/assembly/spades.log",
    conda:
        "../envs/spades.yaml"
    threads: config["assembly_threads"]
    resources:
        mem=config["assembly_memory"],
        time=config["runtime_assembly"],
    shell:
        # remove pipeline_state file to create all output files again
        " rm -f {params.p[outdir]}/pipeline_state/stage_*_copy_files 2> {log} ; "
        " "
        "spades.py "
        " --threads {threads} "
        " --memory {resources.mem} "
        " -o {params.outdir} "
        " -k {params.k}"
        " {params.p[preset]} "
        " {params.p[extra]} "
        " {params.p[inputs]} "
        " {params.p[longreads]} "
        " {params.p[skip_error_correction]} "
        " >> {log} 2>&1 "



# standardizes header labels within contig FASTAs
rule rename_spades_output:
    input:
        rules.run_spades.output
    output:
        multiext("Intermediate/Assembly/spades/{sample}/", "contigs.fasta", "scaffolds.fasta",
        "assembly_graph_with_scaffolds.gfa",
        "scaffolds.paths",
        "contig.paths"),
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