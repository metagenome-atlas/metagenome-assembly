if not config["filter_contigs"]:
    filtered_assembly = raw_assembly

else:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            query=get_quality_controlled_reads,
            target=raw_assembly,
        output:
            bam=temp(
                "Intermediate/Assembly/filtering/{sample}_alignment_to_prefilter_contigs.bam"
            ),
        params:
            extra=config["minimap_extra"],
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log",
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1000,
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

    rule pileup_prefilter:
        input:
            fasta=raw_assembly,
            bam=rules.align_reads_to_prefilter_contigs.output.bam,
        output:
            covstats="Intermediate/Assembly/filtering/{sample}_prefilter_coverage_stats.txt",
        params:
            pileup_secondary="t",
        log:
            "{sample}/logs/assembly/post_process/pilup_prefilter_contigs.log",
        conda:
            "../envs/bbmap.yaml"
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1000,
            java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
        shell:
            "pileup.sh ref={input.fasta} in={input.bam} "
            " threads={threads} "
            " -Xmx{resources.java_mem}G "
            " covstats={output.covstats} "
            " concise=t "
            " secondary={params.pileup_secondary} "
            " 2> {log}"

    rule filter_by_coverage:
        input:
            fasta=raw_assembly,
            covstats="Intermediate/Assembly/filtering/{sample}_prefilter_coverage_stats.txt",
        output:
            fasta=temp("Intermediate/Assembly/filtering/{sample}_final_contigs.fasta"),
        params:
            minc=config["minimum_average_coverage"],
            minp=config["minimum_percent_covered_bases"],
            minr=config["minimum_mapped_reads"],
            minl=config["minimum_contig_length"],
            trim=config["contig_trim_bp"],
        log:
            "{sample}/logs/assembly/post_process/filter_by_coverage.log",
        conda:
            "../envs/bbmap.yaml"
        threads: 1
        resources:
            mem=config["simplejob_mem"],
            java_mem=int(config["simplejob_mem"] * JAVA_MEM_FRACTION),
        shell:
            "filterbycoverage.sh "
            " in={input.fasta} "
            " cov={input.covstats} "
            " out={output.fasta} "
            " minc={params.minc} "
            " minp={params.minp} "
            " minr={params.minr} "
            " minl={params.minl} "
            " trim={params.trim} "
            " -Xmx{resources.java_mem}G "
            " 2> {log}"

    filtered_assembly = rules.filter_by_coverage.output.fasta


# standardizes header labels within contig FASTAs
rule rename_contigs:
    input:
        filtered_assembly,
    output:
        final_assembly,
    conda:
        "../envs/bbmap.yaml"
    threads: config["simplejob_threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime_simplejob"],
    log:
        "logs/assembly/post_process/rename_and_filter_size_{sample}/.log",
    params:
        minlength=config["minimum_contig_length"],
    shell:
        "rename.sh "
        " in={input} out={output} ow=t "
        " prefix={wildcards.sample} "
        " minscaf={params.minlength} &> {log} "


#### contig stats


rule calculate_contigs_stats:
    input:
        lambda wc: final_assembly.format(**wc)
        if (wc.assembly_step == "final")
        else raw_assembly.format(**wc),
    output:
        temp(
            "Intermediate/Assembly/contig_stats/{sample}_{assembly_step}_contig_stats.txt"
        ),
    conda:
        "../envs/bbmap.yaml"
    log:
        "{sample}/logs/assembly/post_process/contig_stats_{assembly_step}.log",
    threads: 1
    resources:
        mem=1,
        time=config["runtime_simplejob"],
    shell:
        "stats.sh in={input} format=3 out={output} &> {log}"


# I dont now if this rule is really necessary
rule combine_sample_contig_stats:
    input:
        expand(
            "Intermediate/Assembly/contig_stats/{{sample}}_{assembly_step}_contig_stats.txt",
            assembly_step=["prefilter", "final"],
        ),
    output:
        "Intermediate/Assembly/contig_stats/{sample}.tsv",
    run:
        import os
        import pandas as pd

        c = pd.DataFrame()
        for f in input:
            df = pd.read_csv(f, sep="\t")
            assembly_step = (
                os.path.basename(f)
                .replace("_contig_stats.txt", "")
                .replace(wildcards.sample + "_", "")
            )
            c.loc[assembly_step]

        c.to_csv(output[0], sep="\t")


### Coverage


rule align_reads_to_final_contigs:
    input:
        query=get_quality_controlled_reads,
        target=final_assembly,
    output:
        bam="Assembly/alignments/{sample}.bam",
    params:
        extra=config["minimap_extra"],
        sorting="coordinate",
    benchmark:
        "logs/benchmarks/assembly/align_reads/{sample}.txt"
    log:
        "logs/assembly/align_reads/{sample}.",
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
    wrapper:
        "v1.19.0/bio/minimap2/aligner"


rule pileup_contigs_sample:
    input:
        fasta=final_assembly,
        bam="Assembly/alignments/{sample}.bam",
    output:
        covhist="Intermediate/Assembly/contig_stats/{sample}_coverage_histogram.txt",
        covstats="Intermediate/Assembly/contig_stats/{sample}_coverage_stats.txt",
        bincov="Intermediate/Assembly/contig_stats/{sample}_coverage_binned.txt",
    params:
        pileup_secondary=("t" if config["count_multi_mapped_reads"] else "f"),
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sample}.txt"
    log:
        "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",  # This log file is uesd for report
    conda:
        "../envs/bbmap.yaml"
    threads: config["threads"]
    resources:
        mem=config["mem"],
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh "
        " ref={input.fasta} "
        " in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " covstats={output.covstats} "
        " hist={output.covhist} "
        " concise=t "
        " secondary={params.pileup_secondary} "
        " bincov={output.bincov} "
        " 2> {log} "


rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    conda:
        "../envs/bbmap.yaml"
    threads: 1
    resources:
        mem=2 * config["simplejob_threads"],
    shell:
        "samtools index {input}"


localrules:
    combine_contig_stats,


rule combine_contig_stats:
    input:
        contig_stats=expand(
            "Intermediate/Assembly/contig_stats/{sample}_final_contig_stats.txt",
            sample=get_all_samples(),
        ),
        gene_tables=expand(
            "Assembly/annotations/genes/{sample}.tsv",
            sample=get_all_samples(),
        ),
        mapping_logs=expand(
            "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",
            sample=get_all_samples(),
        ),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams=expand("Assembly/alignments/{sample}.bam", sample=get_all_samples()),
    output:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    params:
        samples=get_all_samples(),
    log:
        "logs/assembly/combine_contig_stats.log",
    script:
        "../scripts/combine_contig_stats.py"


"""
localrules:
    build_assembly_report,

rule build_assembly_report:
    input:
        combined_contig_stats="stats/combined_contig_stats.tsv",
    output:
        report="reports/assembly_report.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/assembly/report.log",
    script:
        "../report/assembly_report.py"
"""
