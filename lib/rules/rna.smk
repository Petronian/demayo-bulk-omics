rule rna:
    input:
        "results/analysis/overall/feature_counts.txt",
        "results/multiqc/before",
        "results/multiqc/after"

# QC all files before trimming.
rule qc_before:
    input:
        DATA_DIRECTORY + "/{group}"
    output:
        directory("results/qc/before/{group}")
    shell:
        "mkdir -p {output} && fastqc " +
        kwargs2str(combine_kwargs({"-o": "{output}"},
                                  config.get(CONFIG_FASTQC_ARGS, {}))) +
        " {input}/*.fastq*"

# Synthesize results using MultiQC.
rule multiqc_before:
    input:
        expand("results/qc/before/{group}", group=GROUPS)
    output:
        directory("results/multiqc/before")
    shell:
        "mkdir -p {output} && multiqc " +
        kwargs2str(combine_kwargs({"-o": "{output}",
                                   "-f": ""},
                                  config.get(CONFIG_MULTIQC_ARGS, {}))) +
        " {input}"

# QC all files after trimming.
rule qc_after:
    input:
        "results/trimmed_data/{group}/paired" if PAIRED_END else "results/trimmed_data/{group}"
    output:
        directory("results/qc/after/{group}")
    shell:
        "mkdir -p {output} && fastqc " +
        kwargs2str(combine_kwargs({"-o": "{output}"},
                                  config.get(CONFIG_FASTQC_ARGS, {}))) +
        " {input}/*.fastq*"
        
# Synthesize results using MultiQC.
rule multiqc_after:
    input:
        expand("results/qc/after/{group}", group=GROUPS)
    output:
        directory("results/multiqc/after")
    shell:
        "mkdir -p {output} && multiqc " +
        kwargs2str(combine_kwargs({"-o": "{output}",
                                   "-f": ""},
                                  config.get("multiqc-args", {}))) +
        " {input}"

# Binarize and index the SAM files.
rule samtools_namesort:
    input:
        "results/aligned/raw/{file}.sam"
    output:
        temp("results/aligned/tmp/{file}.bam")
    shell:
        "samtools sort " +
        kwargs2str(combine_kwargs({"-o": "{output}",
                                   "-n": ""},
                                  config.get(CONFIG_SAMTOOLS_SORT_ARGS, {}))) +
        " {input}"

# Write mate pair entries needed to deduplicate.
rule samtools_fixmate:
    input:
        "results/aligned/tmp/{file}.bam",
    output:
        temp("results/aligned/tmp/{file}.fixmate.bam"),
    shell:
        "samtools fixmate " +
        kwargs2str(combine_kwargs({"-m": "{input}"},
                                  config.get(CONFIG_SAMTOOLS_FIXMATE_ARGS, {}))) +
        " {output}"

# index the sorted BAM files.
rule samtools_coordsort:
    input:
        "results/aligned/tmp/{file}.fixmate.bam" if DEDUPLICATE else "results/aligned/raw/{file}.sam"
    output:
        temp("results/aligned/tmp/{file}.deduplicate.bam") if DEDUPLICATE else "results/aligned/processed/{file}.bam"
    shell:
        "samtools sort " +
        kwargs2str(combine_kwargs({"-o": "{output}"},
                                  config.get(CONFIG_SAMTOOLS_SORT_ARGS, {}))) +
        " {input}"

# Deduplicate the sorted BAM files.
rule samtools_deduplicate:
    input:
        "results/aligned/tmp/{file}.deduplicate.bam"
    output:
        "results/aligned/processed/{file}.deduplicate.bam"
    shell:
        "samtools markdup " +
        kwargs2str(config.get(CONFIG_SAMTOOLS_MARKDUP_ARGS, {})) +
        " {input} {output}"

# Index the sorted BAM files.
rule samtools_index:
    input:
        "results/aligned/processed/{file}" + input_source + ".bam"
    output:
        "results/aligned/processed/{file}" + input_source + ".bam.csi"
    shell:
        "samtools index " + 
        kwargs2str(combine_kwargs({"-c": "",
                                   "-o": "{output}"},
                                  config.get(CONFIG_SAMTOOLS_MARKDUP_ARGS, {}))) +
        " {input}"

# Trim all files.
rule trim_files:
    input:
        fastq=DATA_DIRECTORY + "/{group}"
    output:
        results=directory("results/trimmed_data/{group}/paired" if PAIRED_END else "results/trimmed_data/{group}"),
        summaryStats="results/trimmed_data/{group}/stats_summary.txt"
    run:
        trim_files(input.fastq,
                   output.results,
                   PAIRED_END,
                   TRIMMOMATIC_TRIMMER, 
                   combine_kwargs({"-summary": output.summaryStats}, config.get(CONFIG_TRIMMOMATIC_ARGS, {})),
                                  threads = JOBLIB_THREADS)

# Align the files with HISAT2.
rule hisat2_align:
    input:
        fastq="results/trimmed_data/{group}/paired" if PAIRED_END else "results/trimmed_data/{group}",
        idx="genome/built_genome" if not ALLOW_PREBUILT_GENOME else ancient("genome/built_genome/")
        # above: allow people to paste in prebuilt genomes
    output:
        results="results/aligned/raw/{group}.sam",
        summaryFile="results/aligned/raw/{group}_summary_file.txt"
    run:
        hisat2_align(input.fastq,
                     input.idx,
                     output.results,
                     PAIRED_END,
                     combine_kwargs({"--summary-file": output.summaryFile}, config.get("hisat2-args", {})))
        
# Generate a count matrix.
rule feature_count:
    input:
        bam=expand("results/aligned/processed/{group}" + input_source + ".bam", group=GROUPS),
        bai=expand("results/aligned/processed/{group}" + input_source + ".bam.csi", group=GROUPS),
        gtf="genome/annotations/" + GENOME_GTF_INFO["rule_fn"]
    output:
        out="results/analysis/overall/feature_counts.txt",
        outSummary="results/analysis/overall/feature_counts.txt.summary"
    shell:
        "featureCounts " + 
        kwargs2str(combine_kwargs({"-a": "{input.gtf}",
                                   "-o": "{output.out}"} | 
                                   ({"-p": "",
                                    "--countReadPairs": ""} if PAIRED_END else {}),
                                  config.get(CONFIG_FEATURECOUNTS_ARGS, {}))) +
        " {input.bam}"