# Bulk -omics data processing

This repository contains information about the bulk processing of -omics data.
This is used by the DeMayo lab
[(RDBL)](https://www.niehs.nih.gov/research/atniehs/labs/rdbl) to process bulk
-omics data.

> [!NOTE]  
> Please open an issue if you notice a bug or have any improvements in mind!
> Any help is appreciated; just don't use the main branch for active development.
> `combined-dev` is the main development branch for this repository; `rnaseq`
> and `chipseq` are for record-keeping and should not be updated further.

## Setting up the pipeline

There are a couple of steps that need to be done. **This tutorial assumes you
are a beginner and uses `mamba` to quickly install all requirements.**

### Download prerequisite programs

1. Ensure `mamba` is installed. You can read instructions
   [here](https://github.com/conda-forge/miniforge) to install it.
   Note that `conda` is equivalent, but slower.

2. Activate `mamba` in your terminal if you haven't already by
   typing `mamba init` into your terminal.

3. Install `git` by typing `mamba install git`.

4. Use `mamba`/`conda` the following to install most necessary prerequisites:
   `mamba install -c bioconda fastqc multiqc trimmomatic hisat2 samtools
   subread deeptools bedtools homer wget snakemake validators`.

5. Use `pip` to install `macs3` while `mamba`/`conda` is activated:
   `pip install macs3`.

6. Install the appropriate genome you need for ChIPseq analysis with HOMER:
   `configureHomer.pl -install [genome name]`.

> [!CAUTION]
> If you install HOMER with `mamba`/`conda`, then `configureHomer.pl` might
> not be added to your command-line `PATH`. You may need to find where the
> file is stored on your system explicitly. For example:
> ```shell
> $ find ~ -name configureHomer.pl
> /ddn/gs1/home/laispf/usr/local/miniforge3/share/homer/configureHomer.pl
> ```
> Then, you can substitute the entire path returned by `find` for
> `configureHomer.pl` in the command above.

After these tasks are completed, you should be able to move onto the next step.

### Set up the `Snakefile` and data to be processed

1. Select a starting point. This can be a new directory or an existing one,
   but do note that all results will be placed at this point.

2. Download this repository such that all files (`Snakefile`, `config.yml`, etc.)
   are placed in your starting folder:
   ```
   your-starting-folder/
      lib/
          ...
      Snakefile
      config.yml
   ```

3. Create a `data` folder:
   ```
   your-starting-folder/
      lib/
          ...
      data/
      Snakefile
      config.yml
   ```

4. Place your files in the `data` folder such that they match the the following
   structure:

    ```
    data/
        group-1/
            file-1.fastq[.gz]
            file-2.fastq[.gz]
            ...
        group-2/
            ...
        control/
            ...
    ...
    ```
    You may have as many experimental groups as you wish, and you may also have
    one control group (recommended name of `control`). The experimental groups
    may have any name. All groups may contain either paired- or single-end
    read data, but this must be consistent across all groups. **Note that the
    control, or input, group is necessary for processing ChIPseq data using
    this pipeline. Also note that multiple files placed within each group
    folder will be treated as multiple technical replicates.**

5. Set up your `config.yml` file (see below). **Ensure that all variables are
   correctly configured for what you hope to analyze.**

6. Execute `snakemake [target] --cores [number of cores]` to begin the
   pipeline. There are three values you can substitute for `[target]`:
   Use `rna` if you simply want a count matrix from your data; use
   `chip` if you want to perform ChIPseq analysis; use `build_genome` if 
   you simply want to build genome files. The number of cores determines
   how many resources on your computer are allocated to this pipeline; I've
   personally noticed that six cores is a good balance between performance
   and resource use.

### Configure your `config.yml` file

Your `config.yml` file stores all of the arguments you wish to pass to the
individual programs in your pipeline to personalize your analysis. Certain
arguments are forbidden from being supplied to the programs since it would
break the pipeline, but many arguments are generally available. The available
arguments are broken down into two key sections:

1. **Pipeline-specific arguments.** These parameters control the behavior of
   the pipeline as a whole and have little to do with the components of the
   pipeline themselves. The table below has a description of the pipeline-
   specific arguments; bold arguments are required/cannot be omitted.

   | Argument | Description |
   | -------- | ----------- |
   | `data-directory` | Location of the data directory. Can be omitted; if so, data directory defaults to `data`. |
   | **`paired-end`** | Are the FASTQ files you're using paired-end? `True` or `False`. |
   | **`genome`** | Genome to use for alignment and feature counting. This can be either (a) the name of a UCSC genome, such as `mm10` or `hg38`, with an optional period-delimited suffix describing the database to take naming from, such as `mm10.knownGene` or `hg38.ensGene`, or (b) paths/links to a FASTQ and GTF file. See the callout below. |
   | `homer-genome` | Genome to use for peak annotation with HOMER. There are two options: (a) if omitted, defaults to the FASTA and GTF files used for the `genome` arguent (corresponding to the files in the `genome/` directory) or (b) if a string, such as `mm10` or `hg38`, corresponds to a dataset installed using HOMER. _If (b), ensure the dataset is installed in HOMER first using `configureHomer.pl` (see above). Custom FASTA and GTF files for HOMER are not supported. |
   | `allow-prebuilt-genome` | If a prebuilt genome (`.ht2l` files) are placed in the `genome/built_genome` folder, do not rebuilt the genome. `True` or (default) `False`. |
   | **`trimmomatic-trimmer`** | Specify the trimmer to use with trimmomatic. Search up online documentation for trimmomatic for details. |
   | `mark-duplicate-reads` | Whether to mark duplicate reads (but not explicitly _remove_ them) using SAMTools. `True` or (default) `False`. |
   | `overall-comparisons` | For ChIPseq, intersect all peak locations and perform annotation, motif analysis, gene ontology analysis, and genome ontology analysis on these intersected peaks. `True` or (default) `False`. |
   | `pairwise-comparisons` | For ChIPseq, perform differential peak calling using `bdgdiff` between all pairs of ChIPseq peak sets. `True` or (default) `False`. |
   | `joblib-threads` | How many threads to allocate to `joblib`-controlled tasks. Any positive integer, default `1`. |

> [!NOTE]
> If providing custom FASTA and GTF files for `genome`, do the following:
>
> ```yaml
> - genome:
>   - fasta: "[link or path to FASTA file]"
>   - gtf: "[link or path to GTF file]"
> ```
>
> Both `fasta` and `gtf` must be specified and not empty strings (`""`). **They can be gzipped (ending with `.gz` or `.gzip`) or uncompressed.**

> [!WARNING]
>
> Custom genomes are not working properly as of right now. Fixes are being implemented.

2. **Program-specific arguments.** These arguments directly control the behavior
   of the programs making up the pipeline. **Please see the individual documentation
   pages of the programs (FastQC, HISAT2, HOMER, etc.) for more information about
   the options available to you.** See the default `config.yml` file provided
   in this repository for an idea of how to use this file.

### Other useful advice

1. **Seeing the commands the pipeline uses:** A good way to see what commands
   the pipeline will execute to process your files is to look at the results of
   a *dry-run*: `snakemake [target] --cores 1 --dry-run -p`. This will
   print the commands that the pipeline will use in **yellow** so you can see 
   what commands/programs are being used as well as the arguments being supplied
   to them.

2. **Seeing a graphical summary of data processing.** This pipeline relies on
   `snakemake`, which offers a graph-based representation of the steps it takes
   to process your files. Ensure that you have `pydot` installed (`mamba install
   pydot`), and then execute `snakemake [rna or chip] --dag | dot -Tsvg > DAG.svg`.
   An SVG file will appear with a graphical summary about the steps the pipeline
   will execute to process your files. See more on
   [`snakemake`'s website.](https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization)
