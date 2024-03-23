## MAKE THESE RULES GENERIC.

# If the genome does not exist, download the hg38 genome.
rule download_genome:
    output:
        "genome/" + GENOME + ".fa"
    shell:
        "mkdir -p genome && wget --continue --tries=0 -O genome/" + GENOME + ".fa.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/" + GENOME + "/bigZips/" + GENOME + ".fa.gz "
        "&& gunzip genome/" + GENOME + ".fa.gz"

# If the annotations do not exist, download the hg38 annotations.
rule download_annotations:
    output:
        "genome/annotations/" + GENOME + ".ncbiRefSeq.gtf"
    shell:
        "mkdir -p genome/annotations && wget --continue --tries=0 -O genome/annotations/" + GENOME + ".ncbiRefSeq.gtf.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/" + GENOME + "/bigZips/genes/" + GENOME + ".ncbiRefSeq.gtf.gz "
        "&& gunzip genome/annotations/" + GENOME + ".ncbiRefSeq.gtf.gz"