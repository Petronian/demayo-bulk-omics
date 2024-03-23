## Snakefile responsible for the entire ChIPseq workflow.
## Combines both the individual processing of groups with the overall analysis later.

# Python imports.
import os
import time
import subprocess
from joblib import Parallel, delayed
from lib.py.validate import find_groups, find_control_group, combine_kwargs, kwargs2str, kwargs2list
from lib.py.strandedness import trim_files, hisat2_align
from lib.py.chip import plot_heatmap, find_intersections, find_all_intersections
from lib.py.names import *

# LOAD IN THE CONFIGURATION FILE.
configfile: "config.yml"

# Create all critical variables, together with defaults.
DATA_DIRECTORY = config.get("data-directory", "data")
GENOME = config["genome"] # critical not to mis-specify so no default
PAIRED_END = config["paired-end"] # critical not to mis-specify so no default
TRIMMOMATIC_TRIMMER = config[CONFIG_TRIMMOMATIC_TRIMMER] # critical not to mis-specify so no default
GROUPS = find_groups(DATA_DIRECTORY)
CONTROL_GROUP = config.get("control-group", find_control_group(GROUPS))
EXPR_GROUPS = [name for name in GROUPS if CONTROL_GROUP is None or name != CONTROL_GROUP]

# Create all additional variables, together with defaults.
tssBaseName = config.get("tss-base-name", "tss")
tssMatrixFileName = "%s_profile_matrix.txt" % tssBaseName
sortedRegionsFileName = "%s_heatmap1_sorted_regions.bed" % tssBaseName
tssGraphFileName = "%s_profile_graph.pdf" % tssBaseName
heatmapGroups = ["All"] + GROUPS
heatmapFilenames = expand("results/analysis/overall/sortUsing{sample}_" + tssGraphFileName, sample = heatmapGroups)

# Rule imports.
include: "lib/rules/chip.smk" # provides "rna" rule
include: "lib/rules/rna.smk" # provides "chip" rule
include: "lib/rules/build.smk"
include: "lib/rules/downloads.smk"
include: "lib/rules/viz.smk"

# Error: you must select a rule to build. Make this PAINFULLY obvious.
rule error:
    run:
        error_str = "You must select a corresponding rule ('rna' or 'chip')."
        for i in range(5): print("!!!!!")
        print("!!!!! NOTE: " + error_str)
        for i in range(5): print("!!!!!")
        raise ValueError("You must select a corresponding rule ('rna' or 'chip').")
