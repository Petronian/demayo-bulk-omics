# Bulk -omics data processing

This repository contains information about the bulk processing of -omics data.
This will eventually be used by the DeMayo lab
[(RDBL)](https://www.niehs.nih.gov/research/atniehs/labs/rdbl) to process bulk
-omics data.

> [!IMPORTANT]  
> No instructions right now while the code is being written. See other
> branches for instructions.

## To-do

You should do these in order; that's why I've ordered things in this way.

- [ ] Separate out rules into `chipseq` and `rnaseq`. Default rule should be
      the more inclusive option, likely RNAseq.
- [ ] Enable the files to process data with multiple replicates. This should
      only involve adapting the existing code for ChIPseq data since the
      existing codebase already does this for RNAseq data. Consider making
      use of the `branch` function for Snakemake.
- [ ] Allow the use of configuration files to set the options available for
      each of the programs that the pipeline uses. Document this.
- [ ] Enable `slurm` support for HPC users so they don't have to wait around.
