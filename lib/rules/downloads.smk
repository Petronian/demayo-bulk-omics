## MAKE THESE RULES GENERIC.

# If the genome does not exist, download the hg38 genome.
rule download_genome:
    output:
        "genome/hg38.fa"
    shell:
        "mkdir -p genome && wget -O genome/hg38.fa.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz "
        "&& gunzip genome/hg38.fa.gz"

# If the annotations do not exist, download the hg38 annotations.
rule download_annotations:
    output:
        "genome/annotations/hg38.ncbiRefSeq.gtf"
    shell:
        "mkdir -p genome/annotations && wget -O genome/annotations/hg38.ncbiRefSeq.gtf.gz "
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz "
        "&& gunzip genome/annotations/hg38.ncbiRefSeq.gtf.gz"