# TODO: add normalization again


rule init_pre_assembly_processing:
    input:  #expect SE or R1,R2 or R1,R2,SE
        get_quality_controlled_reads,
    output:
        temp(
            expand(
                "Intermediate/Assembly/{{sample}}/reads/QC_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input, "in"),
        interleaved="f",
        outputs=lambda wc, output: io_params_for_tadpole(output, "out"),
        verifypaired="t" if PAIRED_END else "f",
    log:
        "{sample}/logs/assembly/init.log",
    conda:
        "../envs/bbmap.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
    shell:
        " reformat.sh "
        " {params.inputs} "
        " interleaved={params.interleaved} "
        " {params.outputs} "
        " iupacToN=t "
        " touppercase=t "
        " qout=33 "
        " overwrite=true "
        " verifypaired={params.verifypaired} "
        " addslash=t "
        " trimreaddescription=t "
        " threads={threads} "
        " pigz=t unpigz=t "
        " -Xmx{resources.java_mem}G "
        " 2> {log} "


rule error_correction:
    input:
        expand(
            "Intermediate/Assembly/{{sample}}/reads/{{previous_steps}}_{fraction}.fastq.gz",
            fraction=MULTIFILE_FRACTIONS,
        ),
    output:
        temp(
            expand(
                "Intermediate/Assembly/{{sample}}/reads/{{previous_steps}}.errorcorr_{fraction}.fastq.gz",
                fraction=MULTIFILE_FRACTIONS,
            )
        ),
    benchmark:
        "logs/benchmarks/assembly/pre_process/{sample}_error_correction_{previous_steps}.txt"
    log:
        "{sample}/logs/assembly/pre_process/error_correction_{previous_steps}.log",
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    params:
        inputs=lambda wc, input: io_params_for_tadpole(input),
        outputs=lambda wc, output: io_params_for_tadpole(output, key="out"),
        prefilter=2,  # Ignore kmers with less than 2 occurance
        minprob=config["error_correction_minprob"],
        tossdepth=config["error_correction_minimum_kmer_depth"],
        tossjunk="t" if config["error_correction_remove_lowdepth"] else "f",
        lowdepthfraction=config["error_correction_lowdepth_fraction"],
        aggressive=config["error_correction_aggressive"],
        shave="f",  # Shave and rinse can produce substantially better assemblies for low-depth data, but they are very slow for large metagenomes.
    threads: config["threads"]
    shell:
        "tadpole.sh -Xmx{resources.java_mem}G "
        " prefilter={params.prefilter} "
        " prealloc=1 "
        " {params.inputs} "
        " {params.outputs} "
        " mode=correct "
        " aggressive={params.aggressive} "
        " tossjunk={params.tossjunk} "
        " lowdepthfraction={params.lowdepthfraction}"
        " tossdepth={params.tossdepth} "
        " merge=t "
        " shave={params.shave} rinse={params.shave} "
        " threads={threads} "
        " pigz=t unpigz=t "
        " ecc=t ecco=t "
        "&> {log} "


rule merge_pairs:
    input:
        unpack(
            lambda wc: {
                fraction: "Intermediate/Assembly/{sample}/reads/{previous_steps}_{fraction}.fastq.gz".format(
                    fraction=fraction, **wc
                )
                for fraction in ["R1", "R2"]
            }
        ),
    output:
        temp(
            expand(
                "Intermediate/Assembly/{{sample}}/reads/{{previous_steps}}.merged_{fraction}.fastq.gz",
                fraction=["R1", "R2", "me"],
            )
        ),
    threads: config["threads"]
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    conda:
        "../envs/bbmap.yaml"
    log:
        "{sample}/logs/assembly/pre_process/merge_pairs_{previous_steps}.log",
    benchmark:
        "logs/benchmarks/assembly/pre_process/merge_pairs_{previous_steps}/{sample}.txt"
    shadow:
        "shallow"
    params:
        kmer=config["merging_k"],
        extend2=config["merging_extend2"],
        flags=config["merging_flags"],
    shell:
        " bbmerge.sh "
        " -Xmx{resources.java_mem}G threads={threads} "
        " in1={input.R1} in2={input.R2} "
        " outmerged={output[2]} "
        " outu={output[0]} outu2={output[1]} "
        " {params.flags} k={params.kmer} "
        " pigz=t unpigz=t "
        " extend2={params.extend2} 2> {log} "
