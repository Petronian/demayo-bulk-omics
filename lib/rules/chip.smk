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
         heatmapFilenames] if OVERALL_COMPARISONS else []

# Generate the TSS matrix using deepTools to generate the TSS heatmap later.
# NOTE that this matrix contains which heatmap values correspond to which sample, so we will only pre-generate the sorted heatmaps
# (rather than sorted matrices) because users can always do this sorting later.
rule tss_matrix:
    input:
        bigWigs=expand("results/coverage/deeptools/{group}.bw", group=GROUPS),
        gtf="genome/annotations/" + GENOME_GTF_INFO["rule_fn"]
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
                                  config.get(CONFIG_COMPUTEMATRIX_REFERENCEPOINT_ARGS, {})))

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
        bam="results/aligned/processed/{group}" + input_source + ".bam",
        bai="results/aligned/processed/{group}" + input_source + ".bam.csi",
        ctrlBam=expand("results/aligned/processed/{group}" + input_source + ".bam", group=CONTROL_GROUP),
        ctrlBai=expand("results/aligned/processed/{group}" + input_source + ".bam.csi", group=CONTROL_GROUP)
    output:
        narrowPeaks="results/analysis/{group}/peaks/macs3_peaks.narrowPeak",
        bed="results/analysis/{group}/peaks/macs3_narrowPeaks.bed",
        bedGraphTreatment="results/analysis/{group}/peaks/macs3_treat_pileup.bdg",
        bedGraphControl="results/analysis/{group}/peaks/macs3_control_lambda.bdg",
        peaksFile="results/analysis/{group}/peaks/macs3_peaks.xls"
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

# Note that this might fail when we have more than two experimental groups.
# rule individual_peak_comparison:
#     input:
#         expand("results/analysis/{group}/peaks/macs3_narrowPeaks.bed", group=EXPR_GROUPS)
#         #expand("results/analysis/{group}/peaks.bed", group=EXPR_GROUPS)
#     output:
#         expand("results/analysis/{group}/peak_intersections.txt", group=EXPR_GROUPS)
#     run:
#         find_intersections(input, output, config.get(CONFIG_BEDTOOLS_INTERSECT_ARGS, {}), threads = JOBLIB_THREADS)

def find_motifs_helper(wildcards):
    base_dict = {"peakFile": f"results/analysis/{wildcards.group}/peaks/macs3_narrowPeaks.bed"}
    if HOMER_CUSTOM: base_dict = base_dict | {"homerGenome": HOMER_GENOME}
    return base_dict

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_motifs:
    input:
        unpack(find_motifs_helper)
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/{group}/motifs")
    shell:
        "findMotifsGenome.pl {input.peakFile} " + HOMER_GENOME + " {output} " +
        kwargs2str(config.get(CONFIG_FINDMOTIFSGENOME_ARGS, {}))

rule annotate_peaks:
    input:
        "results/analysis/{group}/peaks/macs3_narrowPeaks.bed"
        #"results/analysis/{group}/peaks.txt"
    output:
        ann="results/analysis/{group}/annotated_peaks.txt",
        annStats="results/analysis/{group}/annotated_peaks_stats.txt",
        go=directory("results/analysis/{group}/gene_ontology/homer"),
        geno=directory("results/analysis/{group}/genome_ontology/homer"),
    shell:
        "annotatePeaks.pl {input} " + HOMER_GENOME + " " +
        kwargs2str(combine_kwargs({"-go": "{output.go}",
                                   "-genomeOntology": "{output.geno}",
                                   "-annStats": "{output.annStats}"} | 
                                  ({"-gtf": HOMER_ANNOTATION} if HOMER_CUSTOM else {}),
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

def find_combined_motifs_helper(wildcards):
    base_dict = {"peakFile": f"results/analysis/overall/all_intersections.bed"}
    if HOMER_CUSTOM: base_dict = base_dict | {"homerGenome": HOMER_GENOME}
    return base_dict

# Recommended parameter for size: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
rule find_combined_motifs:
    input:
        unpack(find_combined_motifs_helper)
        #"results/analysis/{group}/peaks.txt"
    output:
        directory("results/analysis/overall/motifs_all_intersections")
    shell:
        "findMotifsGenome.pl {input.peakFile} " + HOMER_GENOME + " {output} " +
        kwargs2str(config.get("findMotifsGenome-args", {}))

rule annotate_combined_peaks:
    input:
        "results/analysis/overall/all_intersections.bed"
    output:
        ann="results/analysis/overall/annotated_all_intersections.txt",
        annStats="results/analysis/overall/annotated_all_intersections_stats.txt",
        go=directory("results/analysis/overall/gene_ontology_all_intersections/homer"),
        geno=directory("results/analysis/overall/genome_ontology_all_intersections/homer"),
    shell:
        "annotatePeaks.pl {input} " + HOMER_GENOME + " " + 
        kwargs2str(combine_kwargs({"-go": "{output.go}",
                                   "-genomeOntology": "{output.geno}",
                                   "-annStats": "{output.annStats}"} | 
                                  ({"-gtf": HOMER_ANNOTATION} if HOMER_CUSTOM else {}),
                                  config.get(CONFIG_ANNOTATEPEAKS_ARGS, {}))) +
        " > {output.ann}"