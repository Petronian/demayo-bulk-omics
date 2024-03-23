# Build the HISAT2 index using the hg38 genome.
rule hisat2_build:
    input:
        "genome/" + GENOME + ".fa"
    output:
        directory("genome/built_genome/")
    shell:
        "mkdir -p genome/built_genome/ && hisat2-build " +
        kwargs2str(combine_kwargs({"--large-index": ""},
                                  config.get(CONFIG_HISAT2_BUILD_ARGS, {}))) +
        " {input} genome/built_genome/ref"