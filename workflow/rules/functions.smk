


JAVA_MEM_FRACTION = config["java_mem_fraction"]


def get_preprocessing_steps(config):
    preprocessing_steps = ["QC"]
    if config["normalize_reads_before_assembly"]:
        preprocessing_steps.append("normalized")

    if config["error_correction_before_assembly"]:
        preprocessing_steps.append("errorcorr")

    if config["merge_pairs_before_assembly"] and PAIRED_END:
        preprocessing_steps.append("merged")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)

logger.debug(f"Assembly preprocessing steps: {assembly_preprocessing_steps}")

assembly_params = {}


### general input functions


def io_params_for_tadpole(io, key="in"):
    """This function generates the input flag needed for bbwrap/tadpole for all cases
    possible for get_quality_controlled_reads.

    params:
        io  input or output element from snakemake
        key 'in' or 'out'

        if io contains attributes:
            se -> in={se}
            R1,R2,se -> in1={R1},se in2={R2}
            R1,R2 -> in1={R1} in2={R2}

    """
    N = len(io)
    if N == 1:
        flag = f"{key}1={io[0]}"
    elif N == 2:
        flag = f"{key}1={io[0]} {key}2={io[1]}"
    elif N == 3:
        flag = f"{key}1={io[0]},{io[2]} {key}2={io[1]}"
    else:
        logger.error(
            (
                "File input/output expectation is one of: "
                "1 file = single-end/ interleaved paired-end "
                "2 files = R1,R2, or"
                "3 files = R1,R2,se"
                "got: {n} files:\n{}"
            ).format("\n".join(io), n=len(io))
        )
        sys.exit(1)
    return flag


### workflow functions


def get_all_samples():

    samples= glob_wildcards("{sample}/sequence_quality_control").sample
    return samples