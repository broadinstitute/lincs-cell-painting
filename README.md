# Processed Data for the LINCS Cell Painting Project

The repository stores data and data processing scripts for **a subset** of the [Broad Drug Repurposing Hub](https://clue.io/repurposing#home) collection of compounds.

In this project, the [Connectivity Map](https://clue.io/team) team perturbed A549 cells with ~1,500 compounds across 6 doses in 5 technical replicates. 
We refer to this dataset as `LINCS Pilot 1`.

For a specific list of compounds tested, see [`metadata`](https://github.com/broadinstitute/lincs-cell-painting/tree/master/metadata). 
Information about the compounds can be interactively explored in the [CLUE Repurposing app](https://clue.io/repurposing-app). 
The [Morphology Connectivity Hub](https://clue.io/morphology) is the primary source of this dataset.

## Image-Based Profiling

We apply a unified, image-based profiling pipeline to all 136 384-well plates from `LINCS Pilot 1` .
We use [pycytominer](https://github.com/cytomining/pycytominer) as the primary tool for image-based profiling.

The profiles are processed and stored in the [profiles/](profiles/) directory.
See [`profiles/README.md`](profiles/README.md) for more details and for instructions on how to reproduce the pipeline.

For more details about image-based profiling in general, please refer to [Caicedo et al. 2017](https://doi.org/10.1038/nmeth.4397).

## Computational Environment

We use [conda](https://docs.conda.io/en/latest/) to manage the computational environment.

After installing conda, execute the following to install and navigate to the environment:

```bash
# First, install the `lincs` conda environment
conda env create --force --file environment.yml

# Then, activate the environment and you're all set!
conda activate lincs
```

Also note that when contributing to the repository, make sure to add any new package in the `environment.yml` file.
