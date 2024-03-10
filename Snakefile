## Snakefile responsible for the entire ChIPseq workflow.
## Combines both the individual processing of groups with the overall analysis later.

import os
import time
import subprocess
from joblib import Parallel, delayed

# START USING SETS IN BEDINTERSECT TO INCREASE EFFICIENCY.

# Find all the files containing the groups. At least one should be named
# 'control' (case-insensitive); get rid of that.
GROUPS = [name for name in os.listdir("data/") if os.path.isdir(os.path.join("data", name))]
CONTROL_GROUP = [name for name in GROUPS if name.upper().startswith("CONTROL")][0]
EXPR_GROUPS = [name for name in GROUPS if name != CONTROL_GROUP]

bedGraphBase = "tag_viz"
bedGraphName = bedGraphBase + ".bedGraph"
tssMatrixFileName = "tss_profile_matrix.txt"
sortedRegionsFileName = "tss_heatmap1_sorted_regions.bed"
tssGraphFileName = "tss_profile_graph.pdf"
controlTagDirectory = "results/tag_directores/" + CONTROL_GROUP
heatmapGroups = ["All"] + GROUPS
heatmapFilenames = expand("results/analysis/overall/sortUsing{sample}_" + tssGraphFileName, sample = heatmapGroups)

# Very eloquent of snakemake; we only need motifs and annotated_peaks for the experimental
# groups, but control group is still processed since we need tss_heatmaps for all samples
# including control.
rule all:
    input:
        expand("results/multiqc/before/{group}", group=GROUPS),
        expand("results/multiqc/after/{group}", group=GROUPS),
        expand("results/analysis/{group}/motifs", group=EXPR_GROUPS),
        expand("results/analysis/{group}/annotated_peaks.txt", group=EXPR_GROUPS),
        expand("results/analysis/{group}/peak_intersections.txt", group=EXPR_GROUPS),
        expand("results/analysis/{group}/gene_ontology/homer", group=EXPR_GROUPS),
        expand("results/analysis/{group}/genome_ontology/homer", group=EXPR_GROUPS),
        "results/analysis/overall/annotated_all_intersections.txt",
        "results/analysis/overall/gene_ontology_all_intersections/homer",
        "results/analysis/overall/genome_ontology_all_intersections/homer",
        "results/analysis/overall/motifs_all_intersections",
        heatmapFilenames

# Target: want to have fastqc processing files.
rule counts_all:
    input:
        expand("results/analysis/{group}/feature_counts.txt", group=GROUPS)

# QC all files before trimming.
rule qc_before:
    input:
        "data/{group}"
    output:
        directory("results/qc/before/{group}")
    shell:
        "mkdir -p {output} && fastqc --threads 6 -o {output} {input}/*.fastq*"

# Synthesize results using MultiQC.
rule multiqc_before:
    input:
        "results/qc/before/{group}"
    output:
        directory("results/multiqc/before/{group}")
    shell:
        "mkdir -p {output} && multiqc -f -o {output} {input}"

# Trim all files.
rule trim_files:
    input:
        fastq="data/{group}",
        #qc="results/multiqc/before/{group}"
    output:
        directory("results/trimmed_data/{group}")
    run:
        os.makedirs(str(output))
        argsList = []
        for file in os.listdir(input.fastq):
            argsList.append(["trimmomatic",
                             "SE",
                             "-threads", "6",
                             f"{input.fastq}/{file}",
                             f"{output}/{file}",
                             "SLIDINGWINDOW:51:20"])
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)

# QC all files after trimming.
rule qc_after:
    input:
        "results/trimmed_data/{group}"
    output:
        directory("results/qc/after/{group}")
    shell:
        "mkdir -p {output} && fastqc --threads 6 -o {output} {input}/*.fastq*"

# Synthesize results using MultiQC.
rule multiqc_after:
    input:
        "results/qc/after/{group}"
    output:
        directory("results/multiqc/after/{group}")
    shell:
        "mkdir -p {output} && multiqc -f -o {output} {input}"

# If the genome does not exist, download the mm10 genome.
rule download_genome:
    output:
        "genome/mm10.fa"
    shell:
        "mkdir -p genome && wget -O genome/mm10.fa.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz "
        "&& gunzip genome/mm10.fa.gz"

# If the annotations do not exist, download the mm10 annotations.
rule download_annotations:
    output:
        "genome/annotations/mm10.ncbiRefSeq.gtf"
    shell:
        "mkdir -p genome/annotations && wget -O genome/annotations/mm10.ncbiRefSeq.gtf.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz "
        "&& gunzip genome/annotations/mm10.ncbiRefSeq.gtf.gz"

# Build the HISAT2 index using the mm10 genome.
rule hisat2_build:
    input:
        "genome/mm10.fa"
    output:
        directory("genome/built_genome/")
    shell:
        "mkdir -p genome/built_genome/ && hisat2-build --large-index {input} genome/built_genome/ref"

# Align the files with HISAT2.
rule hisat2_align:
    input:
        fastq="results/trimmed_data/{group}",
        idx="genome/built_genome/",
        #qc="results/qc/after/{group}",
        #multiqc="results/multiqc/after/{group}"
    output:
        "results/aligned/raw/{group}.sam"
    run:
        files = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                          if "fastq" in file])
        subprocess.run(["hisat2",
                        "-p", "6",
                        "--very-sensitive",
                        "--no-spliced-alignment",
                        "-k", "1",
                        "-x", f"{input.idx}/ref",
                        "--no-unal",
                        "-U", f"{files}",
                        "-S", f"{output}"])
    # shell:
    #     "hisat2 --very-sensitive --no-spliced-alignment -k 1 -x {input.idx}/ref --no-unal "
    #     "-U {input.fastq}" + ",".join([file for file in os.listdir(input.fastq)]) + " "
    #     "-S {output}"

# Binarize and index the SAM files.
rule samtools_namesort:
    input:
        "results/aligned/raw/{file}.sam"
    output:
        temp("results/aligned/tmp/{file}.bam")
    shell:
        "samtools sort -@ 6 -n -o {output} {input}"

# Write mate pair entries needed to deduplicate.
rule samtools_fixmate:
    input:
        "results/aligned/tmp/{file}.bam",
    output:
        temp("results/aligned/tmp/{file}.fixmate.bam"),
    shell:
        "samtools fixmate {input} {output}"

# index the sorted BAM files.
rule samtools_coordsort:
    input:
        "results/aligned/tmp/{file}.fixmate.bam"
    output:
        temp("results/aligned/tmp/{file}.deduplicate.bam")
    shell:
        "samtools sort -@ 6 -o {output} {input}"

# Deduplicate the sorted BAM files.
rule samtools_deduplicate:
    input:
        "results/aligned/tmp/{file}.deduplicate.bam"
    output:
        "results/aligned/processed/{file}.deduplicate.bam"
    shell:
        "samtools markdup -@ 6 {input} {output}"

# Index the sorted BAM files.
rule samtools_index:
    input:
        "results/aligned/processed/{file}.deduplicate.bam"
    output:
        "results/aligned/processed/{file}.deduplicate.bam.csi"
    shell:
        "samtools index -c -@ 6 -o {output} {input}"

# Generate a count matrix.
rule feature_count:
    input:
        bam="results/aligned/processed/{group}.deduplicate.bam",
        bai="results/aligned/processed/{group}.deduplicate.bam.csi",
        gtf="genome/annotations/mm10.ncbiRefSeq.gtf"
    output:
        "results/analysis/{group}/feature_counts.txt"
    shell:
        "featureCounts -T 6 -a {input.gtf} -o {output} {input.bam}"

# Create coverage file with extended reads.
rule generate_bam_coverage:
    input:
        bam="results/aligned/processed/{group}.deduplicate.bam",
        bai="results/aligned/processed/{group}.deduplicate.bam.csi"
    output:
        "results/coverage/deeptools/{group}.bw"
    shell:
        "bamCoverage -b {input.bam} -o {output} -of bigwig --effectiveGenomeSize 2652783500 -p 6 --normalizeUsing RPKM -e 300"

# rule bigwig_to_wig:
#     input:
#         "results/coverage/deeptools/{group}.bw"
#     output:
#         "results/coverage/deeptools/{group}.wig"
#     shell:
#         "bigWigToWig {input} {output}"
#     
# # --max-mem or --sort-tmpdir to ease up on memory requirements.
# rule wig_to_bed:
#     input:
#         "results/coverage/deeptools/{group}.wig"
#     output:
#         "results/coverage/deeptools/{group}.bed"
#     shell:
#         "convert2bed -i wig -x < {input} > {output}"

# rule make_tags_control:
#     input:
#         expand("results/aligned/raw/{group}.sam", group=CONTROL_GROUP)
#     output:
#         directory(expand("results/tag_directories/{group}", group=CONTROL_GROUP))
#     shell:
#         "makeTagDirectory {output} -unique -tbp 1 -genome mm10 -checkGC -fragLength 200 {input}"

# Make tag matrices, which also deduplicates.
# Checks for GC bias based on mm10 genome.
# NOTE THAT THIS DOES NOT NORMALIZE. Need to include flag for normalizing if available.
# This should simply only make tag directories for the experimental groups, but
# requires the control directory as a prerequisite. Ensure this occurs in practice.
# rule make_tags_expr:
#     input:
#         alignment="results/aligned/raw/{group}.sam",
#         tag=expand("results/tag_directories/{group}", group=CONTROL_GROUP)
#     output:
#         directory("results/tag_directories/{group}")
#     shell:
#         "makeTagDirectory {output} -unique -tbp 1 -genome mm10 -checkGC -fragLength 200 {input.alignment}"

# Download chrom.sizes file.
# rule download_chrom_sizes:
#     output:
#         "genome/mm10.chrom.sizes"
#     shell:
#         "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"

# Makes a bedGraph file for visualization in UCSC genome browser.
# Note that by not using the -bigWig flag in makeUCSCfile, we won't get track information.
# rule make_ucsc:
#     input:
#         "results/tag_directories/{group}"
#     output:
#         "results/coverage/homer/{group}.bedGraph.gz"
#     shell:
#         "makeUCSCfile {input} -style chipseq -o results/coverage/homer/{wildcards.group}.bedGraph"

# rule bedgraph_to_bigwig:
#     input:
#         gzippedBedGraph="results/coverage/{group}.bedGraph.gz",
#         chromSizes="genome/mm10.chrom.sizes"
#     output:
#         bigWig="results/coverage/{group}.bw",
#         bedGraph="results/coverage/{group}.bedGraph",
#     shell:
#         "gunzip -c {input.gzippedBedGraph} > {output.bedGraph} && bedGraphToBigWig results/coverage/{wildcards.group}.bedGraph "
#         "{input.chromSizes} {output.bigWig}"

# Generate the TSS matrix using deepTools to generate the TSS heatmap later.
# NOTE that this matrix contains which heatmap values correspond to which sample, so we will only pre-generate the sorted heatmaps
# (rather than sorted matrices) because users can always do this sorting later.
rule tss_matrix:
    input:
        bigWigs=expand("results/coverage/deeptools/{group}.bw", group=GROUPS),
        gtf="genome/annotations/mm10.ncbiRefSeq.gtf"
    output:
        mtx="results/analysis/overall/" + tssMatrixFileName,
        gzippedMtx="results/analysis/overall/" + tssMatrixFileName + ".gz",
        sortedRegions="results/analysis/overall/" + sortedRegionsFileName,
    shell:
        "computeMatrix reference-point -S {input.bigWigs} -R {input.gtf} --transcriptID transcript -o {output.gzippedMtx} --outFileNameMatrix {output.mtx} "
        "--outFileSortedRegions {output.sortedRegions} --referencePoint TSS -a 3000 -b 3000 -p 6"

# Generate the TSS heatmap. This is required by the rule 'all'.
# Since I'm making many of these and there's no easy way to parallelize them into arguments, I'm doing it in Python.
rule tss_heatmap:
    input:
        "results/analysis/overall/" + tssMatrixFileName + ".gz"
    output:
        heatmapFilenames
    run:
        argsList = []
        for i, fn in enumerate(heatmapGroups):
            # "All" is first group in the list so this works by making i [1,...,N]
            sortArgs = ["--sortUsingSamples", str(i)] if fn != "All" else []
            argsList.append(["plotHeatmap", "-m", str(input), "-o", str(output[i]), "--colorMap", "Blues"] + sortArgs)
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)
    #shell:
    #    "plotHeatmap -m {input} -o {output}"

## RULES FOR OVERALL RESULTS
# Note that we will need custom Python code from here on out to call shell
# processes with appropriate action for the control.

def find_control(dir):
    dirs_raw = os.listdir(dir)
    for name in os.listdir(dir):
        if name.upper().startswith("CONTROL"):
            return name

# rule find_peaks_homer:
#     input:
#         dir="results/tag_directories/{group}",
#         ucsc="results/coverage/deeptools/{group}.bw",
#         ucsc2="results/coverage/homer/{group}.bedGraph.gz",
#         controlDir=expand("results/tag_directories/{group}", group=CONTROL_GROUP)
#     output:
#         "results/analysis/{group}/peaks.txt"
#     shell:
#         "findPeaks {input.dir} -style factor -o {output} -i {input.controlDir}"

rule find_peaks_macs:
    input:
        #ucsc="results/coverage/deeptools/{group}.bw", # coverage, deeptools; required elsewhere
        #ucsc2="results/coverage/homer/{group}.bedGraph.gz", # coverage, homer; not needed
        bam="results/aligned/processed/{group}.deduplicate.bam",
        bai="results/aligned/processed/{group}.deduplicate.bam.csi",
        ctrlBam=expand("results/aligned/processed/{group}.deduplicate.bam", group=CONTROL_GROUP),
        ctrlBai=expand("results/aligned/processed/{group}.deduplicate.bam.csi", group=CONTROL_GROUP)
    output:
        narrowPeaks="results/analysis/{group}/peaks/macs3_peaks.narrowPeak",
        bed="results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
    shell:
        "macs3 callpeak -t {input.bam} -c {input.ctrlBam} -n macs3 -f BAM "
        "--outdir results/analysis/{wildcards.group}/peaks/ -g mm -q 0.0001 --nomodel --extsize 300 "
        "&& sed -E 's/\..*$/+/g' {output.narrowPeaks} > {output.bed}"

# rule peaks_to_bed:
#     input:
#         "results/analysis/{group}/peaks.txt"
#     output:
#         "results/analysis/{group}/peaks.bed"
#     shell:
#         "pos2bed.pl {input} > {output}"

# Note that this might fail when we have more than two experimental groups.
rule individual_peak_comparison:
    input:
        expand("results/analysis/{group}/peaks/macs3_narrowPeaks.bed", group=EXPR_GROUPS)
        #expand("results/analysis/{group}/peaks.bed", group=EXPR_GROUPS)
    output:
        expand("results/analysis/{group}/peak_intersections.txt", group=EXPR_GROUPS)
    run:
        argList = []
        for i, _ in enumerate(input):
            temp = [str(x) for x in input]
            group = temp.pop(i)
            argList.append((["bedtools", "intersect", "-wo", "-filenames", "-a", group, "-b"] + temp, str(output[i])))
        def find_intersections(argTuple):
            # Write some comments at the top of the file, and then run the bedtools intersect function.
            args, output = argTuple
            file = open(output, "w")
            file.write("# Command: %s\n" % " ".join(args))
            file.write("# See https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html for column information.\n")
            file.close()
            file = open(output, "a")
            subprocess.run(args, text=True, stdout=file)
            file.close()
        Parallel(n_jobs=6)(delayed(find_intersections)(argTuple) for argTuple in argList)

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_motifs:
    input:
        "results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/{group}/motifs")
    shell:
        "findMotifsGenome.pl {input} mm10 {output} -size 200"

rule annotate_peaks:
    input:
        "results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        ann="results/analysis/{group}/annotated_peaks.txt",
        go=directory("results/analysis/{group}/gene_ontology/homer"),
        geno=directory("results/analysis/{group}/genome_ontology/homer"),
    shell:
        "annotatePeaks.pl {input} mm10 -go {output.go} -genomeOntology {output.geno} > {output.ann}"
        
# Note that this might fail when we have more than two experimental groups.
rule overall_peak_comparison:
    input:
        expand("results/analysis/{group}/peaks/macs3_narrowPeaks.bed", group=EXPR_GROUPS)
        #expand("results/analysis/{group}/peaks.bed", group=EXPR_GROUPS)
    output:
        "results/analysis/overall/all_intersections.bed"
    run:
        temp = [str(x) for x in input]
        argsList = ["bedtools", "intersect", "-a", temp[0], "-b"] + temp[1:]
        file = open(str(output), "w")
        subprocess.run(argsList, text=True, stdout=file)
        file.close()

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_combined_motifs:
    input:
        "results/analysis/overall/all_intersections.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/overall/motifs_all_intersections")
    shell:
        "findMotifsGenome.pl {input} mm10 {output} -size 200"

rule annotate_combined_peaks:
    input:
        "results/analysis/overall/all_intersections.bed"
    output:
        ann="results/analysis/overall/annotated_all_intersections.txt",
        go=directory("results/analysis/overall/gene_ontology_all_intersections/homer"),
        geno=directory("results/analysis/overall/genome_ontology_all_intersections/homer"),
    shell:
        "annotatePeaks.pl {input} mm10 -go {output.go} -genomeOntology {output.geno} > {output.ann}"