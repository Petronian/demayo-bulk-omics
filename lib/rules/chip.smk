# Very eloquent of snakemake; we only need motifs and annotated_peaks for the experimental
# groups, but control group is still processed since we need tss_heatmaps for all samples
# including control.
rule chip:
    input:
        expand("results/analysis/{group}/motifs", group=EXPR_GROUPS),
        expand("results/analysis/{group}/annotated_peaks.txt", group=EXPR_GROUPS),
        expand("results/analysis/{group}/peak_intersections.txt", group=EXPR_GROUPS),
        expand("results/analysis/{group}/gene_ontology/homer", group=EXPR_GROUPS),
        expand("results/analysis/{group}/genome_ontology/homer", group=EXPR_GROUPS),
        ["results/analysis/overall/annotated_all_intersections.txt",
         "results/analysis/overall/gene_ontology_all_intersections/homer",
         "results/analysis/overall/genome_ontology_all_intersections/homer",
         "results/analysis/overall/motifs_all_intersections",
         heatmapFilenames] if OVERALL_COMPARISONS else [],
        ["results/analysis/pairwise/differential_peaks"] if PAIRWISE_COMPARISONS else []

# Generate the TSS matrix using deepTools to generate the TSS heatmap later.
# NOTE that this matrix contains which heatmap values correspond to which sample, so we will only pre-generate the sorted heatmaps
# (rather than sorted matrices) because users can always do this sorting later.
rule tss_matrix:
    input:
        bigWigs=expand("results/coverage/deeptools/{group}.bw", group=GROUPS),
        gtf="genome/annotations/" + GENOME + ".ncbiRefSeq.gtf"
    output:
        mtx="results/analysis/overall/" + tssMatrixFileName,
        gzippedMtx="results/analysis/overall/" + tssMatrixFileName + ".gz",
        sortedRegions="results/analysis/overall/" + sortedRegionsFileName,
    shell:
        "computeMatrix reference-point " +
        kwargs2str(combine_kwargs({"-S": "{input.bigWigs}",
                                      "-R": "{input.gtf}",
                                      "--transcriptID": "transcript",
                                      "-o": "{output.gzippedMtx}",
                                      "--outFileNameMatrix": "{output.mtx}",
                                      "--outFileSortedRegions": "{output.sortedRegions}"},
                                  config.get(CONFIG_COMPUTEMATRIX_ARGS, {})))

# Generate the TSS heatmap. This is required by the rule 'all'.
# Since I'm making many of these and there's no easy way to parallelize them into arguments, I'm doing it in Python.
rule tss_heatmap:
    input:
        "results/analysis/overall/" + tssMatrixFileName + ".gz"
    output:
        heatmapFilenames
    run:
        plot_heatmap(input, output, heatmapGroups, config.get(CONFIG_PLOTHEATMAP_ARGS, {}), threads = JOBLIB_THREADS)

## RULES FOR OVERALL RESULTS
# Note that we will need custom Python code from here on out to call shell
# processes with appropriate action for the control.

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
        bed="results/analysis/{group}/peaks/macs3_narrowPeaks.bed",
        bedGraphTreatment="results/analysis/{group}/peaks/macs3_treat_pileup.bdg",
        bedGraphControl="results/analysis/{group}/peaks/macs3_control_lambda.bdg"
    shell:
        "macs3 callpeak " +
        kwargs2str(combine_kwargs({"-t": "{input.bam}",
                                   "-c": "{input.ctrlBam}",
                                   "-f": "BAMPE" if PAIRED_END else "BAM",
                                   "-n": "macs3",
                                   "--outdir": "results/analysis/{wildcards.group}/peaks/",
                                   "--bdg": ""},
                                  config.get(CONFIG_MACS3_CALLPEAK_ARGS, {}))) +
        " && sed -E 's/\..*$/+/g' {output.narrowPeaks} > {output.bed}"

rule pairwise_comparisons_macs:
    input:
        bedGraphTreatment=expand("results/analysis/{group}/peaks/macs3_treat_pileup.bdg", group=EXPR_GROUPS),
        bedGraphControl=expand("results/analysis/{group}/peaks/macs3_control_lambda.bdg", group=EXPR_GROUPS),
        peaksFile=expand("results/analysis/{group}/peaks/macs3_peaks.xls", group=EXPR_GROUPS)
    output:
        directory("results/analysis/pairwise/differential_peaks")
    run:
        pairwise_differential_peakcall(input.bedGraphTreatment,
                                       input.bedGraphControl,
                                       input.peaksFile,
                                       GROUPS,
                                       output,
                                       config.get(CONFIG_MACS3_BDGDIFF_ARGS, {}),
                                       threads = JOBLIB_THREADS)


# Note that this might fail when we have more than two experimental groups.
rule individual_peak_comparison:
    input:
        expand("results/analysis/{group}/peaks/macs3_narrowPeaks.bed", group=EXPR_GROUPS)
        #expand("results/analysis/{group}/peaks.bed", group=EXPR_GROUPS)
    output:
        expand("results/analysis/{group}/peak_intersections.txt", group=EXPR_GROUPS)
    run:
        find_intersections(input, output, config.get(CONFIG_BEDTOOLS_INTERSECT_ARGS, {}), threads = JOBLIB_THREADS)

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_motifs:
    input:
        "results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/{group}/motifs")
    shell:
        "findMotifsGenome.pl {input} " + GENOME + " {output} " +
        kwargs2str(config.get(CONFIG_FINDMOTIFSGENOME_ARGS, {}))

rule annotate_peaks:
    input:
        "results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        ann="results/analysis/{group}/annotated_peaks.txt",
        go=directory("results/analysis/{group}/gene_ontology/homer"),
        geno=directory("results/analysis/{group}/genome_ontology/homer"),
    shell:
        "annotatePeaks.pl {input} " + GENOME + " " +
        kwargs2str(combine_kwargs({"-go": "{output.go}",
                                   "-genomeOntology": "{output.geno}"},
                                  config.get(CONFIG_ANNOTATEPEAKS_ARGS, {}))) +
        " > {output.ann}"
        
# Note that this might fail when we have more than two experimental groups.
rule overall_peak_comparison:
    input:
        expand("results/analysis/{group}/peaks/macs3_narrowPeaks.bed", group=EXPR_GROUPS)
        #expand("results/analysis/{group}/peaks.bed", group=EXPR_GROUPS)
    output:
        "results/analysis/overall/all_intersections.bed"
    run:
        find_all_intersections(input, output, config.get(CONFIG_BEDTOOLS_INTERSECT_ARGS, {}))

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_combined_motifs:
    input:
        "results/analysis/overall/all_intersections.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/overall/motifs_all_intersections")
    shell:
        "findMotifsGenome.pl {input} " + GENOME + " {output} " +
        kwargs2str(config.get("findMotifsGenome-args", {}))

rule annotate_combined_peaks:
    input:
        "results/analysis/overall/all_intersections.bed"
    output:
        ann="results/analysis/overall/annotated_all_intersections.txt",
        go=directory("results/analysis/overall/gene_ontology_all_intersections/homer"),
        geno=directory("results/analysis/overall/genome_ontology_all_intersections/homer"),
    shell:
        "annotatePeaks.pl {input} " + GENOME + " " + 
        kwargs2str(combine_kwargs({"-go": "{output.go}",
                                   "-genomeOntology": "{output.geno}"},
                                  config.get(CONFIG_ANNOTATEPEAKS_ARGS, {}))) +
        " > {output.ann}"