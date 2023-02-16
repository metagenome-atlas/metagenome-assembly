# standardizes header labels within contig FASTAs
rule rename_contigs:
    input:
        raw_assembly,
    output:
        "Intermediate/Assembly/{sample}/{sample}_prefilter_contigs.fasta",
    conda:
        "../envs/bbmap.yaml"
    threads: config.get("simplejob_threads", 1)
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime_simplejob"],
    log:
        "{sample}/logs/assembly/post_process/rename_and_filter_size.log",
    params:
        minlength=config["minimum_contig_length"],
    shell:
        "rename.sh "
        " in={input} out={output} ow=t "
        " prefix={wildcards.sample} "
        " minscaf={params.minlength} &> {log} "


rule calculate_contigs_stats:
    input:
        "Intermediate/Assembly/{sample}/{sample}_{assembly_step}_contigs.fasta",
    output:
        "Intermediate/Assembly/{sample}/contig_stats/{assembly_step}_contig_stats.txt",
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


rule combine_sample_contig_stats:
    input:
        expand(
            "Intermediate/Assembly/{{sample}}/contig_stats/{assembly_step}_contig_stats.txt",
            assembly_step=["prefilter", "final"],
        ),
    output:
        "Intermediate/Assembly/{sample}/contig_stats.tsv",
    run:
        import os
        import pandas as pd

        c = pd.DataFrame()
        for f in input:
            df = pd.read_csv(f, sep="\t")
            assembly_step = os.path.basename(f).replace("_contig_stats.txt", "")
            c.loc[assembly_step]

        c.to_csv(output[0], sep="\t")


if config["filter_contigs"]:

    ruleorder: align_reads_to_prefilter_contigs > align_reads_to_final_contigs

    rule align_reads_to_prefilter_contigs:
        input:
            query=get_quality_controlled_reads,
            target=rules.rename_contigs.output,
        output:
            bam=temp("{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam"),
        params:
            extra="-x sr",
        log:
            "{sample}/logs/assembly/post_process/align_reads_to_prefiltered_contigs.log",
        threads: config["threads"]
        resources:
            mem_mb=config["mem"] * 1000,
        wrapper:
            "v1.19.0/bio/minimap2/aligner"

    rule pileup_prefilter:
        input:
            fasta="Intermediate/Assembly/{sample}/{sample}_prefilter_contigs.fasta",
            bam="{sample}/sequence_alignment/alignment_to_prefilter_contigs.bam",
        output:
            covstats="Intermediate/Assembly/{sample}/contig_stats/prefilter_coverage_stats.txt",
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
            fasta="Intermediate/Assembly/{sample}/{sample}_prefilter_contigs.fasta",
            covstats="Intermediate/Assembly/{sample}/contig_stats/prefilter_coverage_stats.txt",
        output:
            fasta="Intermediate/Assembly/{sample}/{sample}_final_contigs.fasta",
            removed_names="Intermediate/Assembly/{sample}/{sample}_discarded_contigs.fasta",
        params:
            minc=config["minimum_average_coverage"],
            minp=config["minimum_percent_covered_bases"],
            minr=config.get("minimum_mapped_reads", MINIMUM_MAPPED_READS),
            minl=config.get("minimum_contig_length", MINIMUM_CONTIG_LENGTH),
            trim=config.get("contig_trim_bp", CONTIG_TRIM_BP),
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
            " outd={output.removed_names} "
            " minc={params.minc} "
            " minp={params.minp} "
            " minr={params.minr} "
            " minl={params.minl} "
            " trim={params.trim} "
            " -Xmx{resources.java_mem}G "
            " 2> {log}"


# HACK: this makes two copies of the same file


else:  # no filter

    localrules:
        do_not_filter_contigs,

    rule do_not_filter_contigs:
        input:
            rules.rename_contigs.output,
        output:
            "Intermediate/Assembly/{sample}/{sample}_final_contigs.fasta",
        threads: 1
        shell:
            "cp {input} {output}"


localrules:
    finalize_contigs,


rule finalize_contigs:
    input:
        "Intermediate/Assembly/{sample}/{sample}_final_contigs.fasta",
    output:
        "{sample}/{sample}_contigs.fasta",
    threads: 1
    run:
        os.symlink(os.path.relpath(input[0], os.path.dirname(output[0])), output[0])


# generalized rule so that reads from any "sample" can be aligned to contigs from "sample_contigs"
rule align_reads_to_final_contigs:
    input:
        query=get_quality_controlled_reads,
        target="{sample_contigs}/{sample_contigs}_contigs.fasta",
    output:
        bam="{sample_contigs}/sequence_alignment/{sample}.bam",
    params:
        extra="-x sr",
        sorting="coordinate",
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/align_reads_to_filtered_contigs/{sample}_to_{sample_contigs}.txt"
    log:
        "{sample_contigs}/logs/assembly/calculate_coverage/align_reads_from_{sample}_to_filtered_contigs.log",
    threads: config["threads"]
    resources:
        mem_mb=config["mem"] * 1000,
    wrapper:
        "v1.19.0/bio/minimap2/aligner"


rule pileup_contigs_sample:
    input:
        fasta="{sample}/{sample}_contigs.fasta",
        bam="{sample}/sequence_alignment/{sample}.bam",
    output:
        covhist="Intermediate/Assembly/{sample}/contig_stats/postfilter_coverage_histogram.txt",
        covstats="Intermediate/Assembly/{sample}/contig_stats/postfilter_coverage_stats.txt",
        bincov="Intermediate/Assembly/{sample}/contig_stats/postfilter_coverage_binned.txt",
    params:
        pileup_secondary=(
            "t"
            if config["count_multi_mapped_reads"]
            else "f"
        ),
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sample}.txt"
    log:
        "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",  # This log file is uesd for report
    conda:
        "../envs/bbmap.yaml"
    threads: config.get("threads", 1)
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


rule predict_genes:
    input:
        "{sample}/{sample}_contigs.fasta",
    output:
        fna="{sample}/annotation/predicted_genes/{sample}.fna",
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
        gff="{sample}/annotation/predicted_genes/{sample}.gff",
    conda:
        "../envs/prodigal.yaml"
    log:
        "{sample}/logs/gene_annotation/prodigal.txt",
    benchmark:
        "logs/benchmarks/prodigal/{sample}.txt"
    threads: 1
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime_simplejob"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """


localrules:
    get_contigs_from_gene_names,


rule get_contigs_from_gene_names:
    input:
        faa="{sample}/annotation/predicted_genes/{sample}.faa",
    output:
        tsv="{sample}/annotation/predicted_genes/{sample}.tsv",
    run:
        header = [
            "gene_id",
            "Contig",
            "Gene_nr",
            "Start",
            "Stop",
            "Strand",
            "Annotation",
        ]
        with open(output.tsv, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            with open(input.faa) as fin:
                gene_idx = 0
                for line in fin:
                    if line[0] == ">":
                        text = line[1:].strip().split(" # ")
                        old_gene_name = text[0]
                        text.remove(old_gene_name)
                        old_gene_name_split = old_gene_name.split("_")
                        gene_nr = old_gene_name_split[-1]
                        contig_nr = old_gene_name_split[-2]
                        sample = "_".join(
                            old_gene_name_split[: len(old_gene_name_split) - 2]
                        )
                        tsv.write(
                            "{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                                text="\t".join(text),
                                gene_id=old_gene_name,
                                i=gene_idx,
                                sample=sample,
                                gene_nr=gene_nr,
                                contig_nr=contig_nr,
                            )
                        )
                        gene_idx += 1
                        #



localrules:
    combine_contig_stats,


rule combine_contig_stats:
    input:
        contig_stats=expand(
            "Intermediate/Assembly/{sample}/contig_stats/final_contig_stats.txt", sample=get_all_samples()
        ),
        gene_tables=expand(
            "{sample}/annotation/predicted_genes/{sample}.tsv", sample=get_all_samples()
        ),
        mapping_logs=expand(
            "{sample}/logs/assembly/calculate_coverage/pilup_final_contigs.log",
            sample=get_all_samples(),
        ),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams=expand("{sample}/sequence_alignment/{sample}.bam", sample=get_all_samples()),
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
