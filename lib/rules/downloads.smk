## MAKE THESE RULES GENERIC.

# If the genome does not exist, download the hg38 genome.
rule download_genome:
    output:
        "genome/" + GENOME_FASTA_RULE_NAME
    shell:
        "mkdir -p genome/ && " +
        ("wget --continue --tries 0 -O genome/" + GENOME_FASTA_INFO["fn"] +
         " " + GENOME_FASTA_INFO["fp"] if GENOME_FASTA_INFO["is_url"]
                                       else "cp " + GENOME_FASTA_INFO["fp"] + " genome/") +
        (" && gunzip genome/" + GENOME_FASTA_INFO["fn"] if GENOME_FASTA_INFO["gzipped"]
                                                        else "")

# If the annotations do not exist, download the hg38 annotations.
rule download_annotations:
    output:
        "genome/annotations/" + GENOME_GTF_RULE_NAME
    shell:
        "mkdir -p genome/annotations/ && " +
        ("wget --continue --tries 0 -O genome/annotations/" + GENOME_GTF_INFO["fn"] +
         " " + GENOME_GTF_INFO["fp"] if GENOME_GTF_INFO["is_url"]
                                     else "cp " + GENOME_GTF_INFO["fp"] + " genome/annotations/") +
        (" && gunzip genome/annotations/" + GENOME_GTF_INFO["fn"] if GENOME_GTF_INFO["gzipped"]
                                                                 else "")