# Build the HISAT2 index.
rule build_genome:
    input:
        "genome/" + GENOME_FASTA_INFO["rule_fn"] if not ALLOW_PREBUILT_GENOME
                                                 else ancient("genome/" + GENOME_FASTA_INFO["rule_fn"])
        # above: allow people to paste in prebuilt genomes
    output:
        directory("genome/built_genome/")
    shell:
        "mkdir -p genome/built_genome/ && hisat2-build " +
        kwargs2str(combine_kwargs({"--large-index": ""},
                                  config.get(CONFIG_HISAT2_BUILD_ARGS, {}))) +
        " {input} genome/built_genome/ref"