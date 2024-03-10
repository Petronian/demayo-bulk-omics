# Bulk RNAseq data processing

This folder contains information about the bulk processing of RNAseq data for the GD 2.5 oviduct project for SRF/MYOCD mice.

## Step 1: create your data

Navigate to the directory containing the `Snakefile` and create a folder named `data` (if not already present). This pipeline takes in FASTQ files containing raw sequencing data. The data need to be placed in the appropriate directory structure before proceeding:

```
data/
    group-1/
         file-1.fastq[.gz]
         file-2.fastq[.gz]
         ...
    group-2/
         ...
    ...
```

Each of the FASTQ files in the group directories will be treated as ONE REPLICATE, simply concatenating all of the relevant reads to an associated SAM file (which will then be processed).

## Step 2: run `snakemake`

Again navigate to the directory containing the `Snakefile` and execute the following command:

```sh
snakemake -p --cores [cores]
```

This will run the analysis using `[cores]` number of CPUs. Let this go, and eventually you should have your results.