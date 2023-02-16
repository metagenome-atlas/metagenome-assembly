'''

PAIRED_END = sampleTable.columns.str.contains("R2").any() or config.get(
    "interleaved_fastqs", False
)


colum_headers_QC = sampleTable.columns[sampleTable.columns.str.startswith("Reads_QC_")]
if len(colum_headers_QC) >= 1:
    MULTIFILE_FRACTIONS = list(colum_headers_QC.str.replace("Reads_QC_", ""))

    if (len(MULTIFILE_FRACTIONS) == 1) and config.get("interleaved_fastqs", False):
        MULTIFILE_FRACTIONS = ["R1", "R2"]

else:
    MULTIFILE_FRACTIONS = ["R1", "R2", "se"] if PAIRED_END else ["se"]

colum_headers_raw = sampleTable.columns[
    sampleTable.columns.str.startswith("Reads_raw_")
]
if len(colum_headers_raw) == 0:
    SKIP_QC = True

    logger.info("Didn't find raw reads in sampleTable - skip QC")
    RAW_INPUT_FRACTIONS = MULTIFILE_FRACTIONS
else:
    RAW_INPUT_FRACTIONS = ["R1", "R2"] if PAIRED_END else ["se"]


if (len(colum_headers_raw) == 0) and (len(colum_headers_QC) == 0):
    raise IOError(
        "Either raw reas or QC reads need to be in the sample table. "
        "I din't find any columnns with 'Reads_raw_<fraction>' or 'Reads_QC_<fraction>'  "
    )


def get_quality_controlled_reads(wildcards, include_se=False):
    """
    Gets quality controlled reads.
    R1 and R1 or se are returned as a dict.

    if the files are not in the sample tible impute default path produced with atlas.
    set

    """

    Fractions = MULTIFILE_FRACTIONS

    if config.get("interleaved_fastqs", False) and SKIP_QC:
        Fractions = ["se"]

    elif not include_se:
        # get only R1 and R2 or se
        Fractions = Fractions[: min(len(Fractions), 2)]

    try:
        QC_Headers = ["Reads_QC_" + f for f in Fractions]
        return get_files_from_sampleTable(wildcards.sample, QC_Headers)

    except FileNotInSampleTableException:
        # return files as named by atlas pipeline
        return expand(
            "{sample}/sequence_quality_control/{sample}_QC_{fraction}.fastq.gz",
            fraction=Fractions,
            sample=wildcards.sample,
        )

'''
## new

PAIRED_END = True
MULTIFILE_FRACTIONS= ["R1", "R2", "se"] if PAIRED_END else ["se"]
JAVA_MEM_FRACTION = 0.85

def get_quality_controlled_reads(wildcards):

    return expand(
        "{sample}/sequence_quality_control/{sample}_QC_{fraction}.fastq.gz",
        fraction=MULTIFILE_FRACTIONS,
        sample=wildcards.sample,
    )


