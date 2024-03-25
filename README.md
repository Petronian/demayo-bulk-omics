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
   subread deeptools bedtools homer wget snakemake`.

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

6. Execute `snakemake [rna or chip] --cores [number of cores]` to begin the
   pipeline. Use `rna` if you simply want a count matrix from your data; use
   `chip` if you want to perform ChIPseq analysis. The number of cores determines
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
   pipeline themselves.

2. **Program-specific arguments.** These arguments directly control the behavior
   of the programs making up the pipeline. **Please see the individual documentation
   pages of the programs (FastQC, HISAT2, HOMER, etc.) for more information about
   the options available to you.** See the default `config.yml` file provided
   in this repository for an idea of how to use this file.

### Other useful advice

1. **Seeing the commands the pipeline uses:** A good way to see what commands
   the pipeline will execute to process your files is to look at the results of
   a *dry-run*: `snakemake [rna or chip] --cores 1 --dry-run -p`. This will
   print the commands that the pipeline will use in **yellow** so you can see 
   what commands/programs are being used as well as the arguments being supplied
   to them.

2. **Seeing a graphical summary of data processing.** This pipeline relies on
   `snakemake`, which offers a graph-based representation of the steps it takes
   to process your files. Ensure that you have `pydot` installed (`mamba install
   pydot`), and then execute `snakemake [rna or chip] --dag | dot -Tsvg > DAG.svg`.
   An SVG file will appear with a graphical summary about the steps the pipeline
   will execute to process your files.See more on
   [`snakemake`'s website.](https://snakemake.readthedocs.io/en/stable/executing/cli.html#visualization)

## To-do

These are ordered tasks that I wish to complete before finalizing the pipeline.

- [x] Separate out rules into `chipseq` and `rnaseq`. Default rule should be
      the more inclusive option, likely RNAseq.
- [x] Enable the files to process data with multiple replicates. This should
      only involve adapting the existing code for ChIPseq data since the
      existing codebase already does this for RNAseq data. Consider making
      use of the `branch` function for Snakemake.
- [x] Allow the use of configuration files to set the options available for
      each of the programs that the pipeline uses. Document this.
- [ ] Enable `slurm` support for HPC users so they don't have to wait around.
- [ ] Clean up the atrocious coding patterns used in this repository.
