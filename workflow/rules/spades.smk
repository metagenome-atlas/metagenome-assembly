
from copy import deepcopy

if PAIRED_END:
    ASSEMBLY_FRACTIONS = ["R1", "R2"]
    if config.get("merge_pairs_before_assembly", True):
        ASSEMBLY_FRACTIONS += ["me"]
else:
    ASSEMBLY_FRACTIONS = deepcopy(MULTIFILE_FRACTIONS)

    if config["spades_preset"] == "meta":
        logger.error(
            "Metaspades cannot handle single end libraries. Use another assembler or specify 'spades_preset': normal"
        )
        exit(1)

assembly_params["spades"] = {"meta": "--meta", "normal": "", "rna": "--rna"}


def spades_parameters(wc, input):
    if not os.path.exists("Intermediate/Assembly/{sample}/params.txt".format(sample=wc.sample)):
        params = {}

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
            long_read_file = get_files_from_sampleTable(wc.sample, "longreads")[0]
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

    else:
        params = {
            "inputs": "--restart-from last",
            "preset": "",
            "skip_error_correction": "",
            "extra": "",
            "longreads": "",
        }

    params["outdir"] = "Intermediate/Assembly/{sample}/spades".format(sample=wc.sample)

    return params


rule run_spades:
    input:
        expand(
            "Intermediate/Assembly/{{sample}}/reads/{assembly_preprocessing_steps}_{fraction}.fastq.gz",
            fraction=ASSEMBLY_FRACTIONS,
            assembly_preprocessing_steps=assembly_preprocessing_steps,
        ),
    output:
        "Intermediate/Assembly/{sample}/spades/contigs.fasta",
        "Intermediate/Assembly/{sample}/spades/scaffolds.fasta",
    benchmark:
        "logs/benchmarks/assembly/spades/{sample}.txt"
    params:
        p=lambda wc, input: spades_parameters(wc, input),
        k=config["spades_k"]
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
        " -o {params.p[outdir]} "
        " -k {params.k}"
        " {params.p[preset]} "
        " {params.p[extra]} "
        " {params.p[inputs]} "
        " {params.p[longreads]} "
        " {params.p[skip_error_correction]} "
        " >> {log} 2>&1 "
