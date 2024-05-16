## Snakefile responsible for the entire ChIPseq workflow.
## Combines both the individual processing of groups with the overall analysis later.

# Python imports.
import os
import time
import subprocess
from joblib import Parallel, delayed
from lib.py.validate import find_groups, find_control_group, combine_kwargs, kwargs2str, kwargs2list, process_genome_argument, process_homer_genome_argument
from lib.py.strandedness import trim_files, hisat2_align
from lib.py.chip import plot_heatmap, find_intersections, find_all_intersections, pairwise_differential_peakcall
from lib.py.names import *

# LOAD IN THE CONFIGURATION FILE.
configfile: "config.yml"

# Create all critical variables, together with defaults.
DATA_DIRECTORY = config.get("data-directory", "data")
GENOME = config["genome"] # critical not to mis-specify so no default
ALLOW_PREBUILT_GENOME = config.get("allow-prebuilt-genome", False) # allow pasting of prebuilt genome
PAIRED_END = config["paired-end"] # critical not to mis-specify so no default
OVERALL_COMPARISONS = config.get("overall-comparisons", False) # critical not to mis-specify so no default
TRIMMOMATIC_TRIMMER = config[CONFIG_TRIMMOMATIC_TRIMMER] # critical not to mis-specify so no default
GROUPS = find_groups(DATA_DIRECTORY)
CONTROL_GROUP = config.get("control-group", find_control_group(GROUPS))
EXPR_GROUPS = [name for name in GROUPS if CONTROL_GROUP is None or name != CONTROL_GROUP]
JOBLIB_THREADS = int(config.get("joblib-threads", 1))
DEDUPLICATE = config.get("mark-duplicate-reads", False)

# A special spot to process genomes to be used.
GENOME_FASTA_INFO, GENOME_GTF_INFO = process_genome_argument(config["genome"])
GENOME_FASTA_RULE_NAME = GENOME_FASTA_INFO["rule_fn"]
GENOME_GTF_RULE_NAME = GENOME_GTF_INFO["rule_fn"]
HOMER_GENOME, HOMER_ANNOTATION, HOMER_CUSTOM = process_homer_genome_argument(
    config.get("homer-genome"), GENOME_FASTA_INFO, GENOME_GTF_INFO)

# Create all additional variables, together with defaults.
tssBaseName = config.get("tss-base-name", "tss")
tssMatrixFileName = "%s_profile_matrix.txt" % tssBaseName
sortedRegionsFileName = "%s_heatmap1_sorted_regions.bed" % tssBaseName
tssGraphFileName = "%s_profile_graph.pdf" % tssBaseName
heatmapGroups = ["All"] + GROUPS
heatmapFilenames = expand("results/analysis/overall/sortUsing{sample}_" + tssGraphFileName, sample = heatmapGroups)
input_source = ".deduplicate" if DEDUPLICATE else "" # Remove duplicates or don't

# Rule imports.
include: "lib/rules/chip.smk" # provides "rna" rule
include: "lib/rules/rna.smk" # provides "chip" rule
include: "lib/rules/build.smk" # provides "build_genome" rule
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
