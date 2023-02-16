import os
import re
import sys
from glob import glob
from snakemake.utils import report
import warnings
from copy import deepcopy


def get_preprocessing_steps(config):
    preprocessing_steps = ["QC"]
    if config.get("normalize_reads_before_assembly", False):
        preprocessing_steps.append("normalized")

    if config.get("error_correction_before_assembly", True):
        preprocessing_steps.append("errorcorr")

    if config.get("merge_pairs_before_assembly", True) and PAIRED_END:
        preprocessing_steps.append("merged")

    return ".".join(preprocessing_steps)


assembly_preprocessing_steps = get_preprocessing_steps(config)


assembly_params = {}