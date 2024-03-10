## Snakefile responsible for the entire ChIPseq workflow.
## Combines both the individual processing of groups with the overall analysis later.

import os
import time
import subprocess
from joblib import Parallel, delayed

GROUPS = [name for name in os.listdir("data/") if os.path.isdir(os.path.join("data", name))]

# Target: want to have fastqc processing files.
rule counts_all:
    input:
        "results/analysis/overall/feature_counts.txt",
        "results/multiqc/before",
        "results/multiqc/after"

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
        expand("results/qc/before/{group}", group=GROUPS)
    output:
        directory("results/multiqc/before")
    shell:
        "mkdir -p {output} && multiqc -f -o {output} {input}"

# Trim all files.
rule trim_files:
    input:
        fastq="data/{group}",
        #qc="results/multiqc/before/{group}"
    output:
        directory("results/trimmed_data/{group}/paired")
    run:
        # Preparatory work.
        parent = os.path.abspath(os.path.join(str(output), os.pardir))
        os.makedirs(parent + "/paired", exist_ok = True)
        os.makedirs(parent + "/unpaired", exist_ok = True)
        argsList = []
        
        # Sort files into 1's and 2's.
        file1s = {}
        file2s = {}
        for file in os.listdir(input.fastq):
            if ".1.fastq" in file:
                substr_start = file.find(".1.fastq")
                file1s[file[0:(substr_start - 1)]] = file
            elif ".2.fastq" in file:
                substr_start = file.find(".2.fastq")
                file2s[file[0:(substr_start - 1)]] = file
            else:
                raise ValueError("File %s is not part of a pair." % file)
        
        # Now go through all the files we have.
        basenames = set(file1s.keys()) | set(file2s.keys())
        for basename in basenames:
            try:
                file1 = file1s[basename]
                file2 = file2s[basename]
            except ValueError:
                raise ValueError("There must be a file that doesn't have "
                                 "both paired-end reads (.1.fastq and .2.fastq). "
                                 "Check your files.")
            argsList.append(["trimmomatic",
                             "PE",
                             "-threads", "6",
                             f"{input.fastq}/{file1}",
                             f"{input.fastq}/{file2}",
                             f"{parent}/paired/{basename}.paired.1.fastq.gz",
                             f"{parent}/unpaired/{basename}.unpaired.1.fastq.gz",
                             f"{parent}/paired/{basename}.paired.2.fastq.gz",
                             f"{parent}/unpaired/{basename}.unpaired.2.fastq.gz",
                             "SLIDINGWINDOW:51:20"])
        Parallel(n_jobs=6)(delayed(subprocess.run)(jobArgSet) for jobArgSet in argsList)

# QC all files after trimming.
rule qc_after:
    input:
        "results/trimmed_data/{group}/paired"
    output:
        directory("results/qc/after/{group}")
    shell:
        "mkdir -p {output} && fastqc --threads 6 -o {output} {input}/*.fastq*"

# Synthesize results using MultiQC.
rule multiqc_after:
    input:
        expand("results/qc/after/{group}", group=GROUPS)
    output:
        directory("results/multiqc/after")
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
rule hisat2_align_paired_end:
    input:
        fastq="results/trimmed_data/{group}/paired",
        idx="genome/built_genome/",
        #qc="results/qc/after/{group}",
        #multiqc="results/multiqc/after/{group}"
    output:
        "results/aligned/raw/{group}.sam"
    run:
        files_1 = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                           if ".1.fastq" in file])
        files_2 = ",".join([os.path.join(input.fastq, file) for file in os.listdir(input.fastq)
                           if ".2.fastq" in file])
        subprocess.run(["hisat2",
                        "-p", "6",
                        "--very-sensitive",
                        "-k", "1",
                        "-x", f"{input.idx}/ref",
                        "--no-unal",
                        "-1", f"{files_1}",
                        "-2", f"{files_2}",
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
        "samtools fixmate -m {input} {output}"

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
        bam=expand("results/aligned/processed/{group}.deduplicate.bam", group=GROUPS),
        bai=expand("results/aligned/processed/{group}.deduplicate.bam.csi", group=GROUPS),
        gtf="genome/annotations/mm10.ncbiRefSeq.gtf"
    output:
        "results/analysis/overall/feature_counts.txt"
    shell:
        "featureCounts -T 6 -a {input.gtf} -o {output} {input.bam}"
