# Create coverage file with extended reads.
rule generate_bam_coverage:
    input:
        bam="results/aligned/processed/{group}.deduplicate.bam",
        bai="results/aligned/processed/{group}.deduplicate.bam.csi"
    output:
        "results/coverage/deeptools/{group}.bw"
    shell:
        "bamCoverage " +
        kwargs2str(combine_kwargs({"-b": "{input.bam}",
                                   "-o": "{output}",
                                   "-of": "bigwig"},
                                  config.get(CONFIG_BAMCOVERAGE_ARGS, {})))

