# Processed Cell Painting Data for the LINCS Drug Repurposing Project

The repository stores data and data processing scripts for **a subset** of the LINCS Drug Repurposing Project.
In the LINCS Cell Painting project, we perturbed A549 cells with ~1,500 compounds across 6 doses with ~4 biological replicates.

For a specific list of compounds tested, see [`metadata`](https://github.com/broadinstitute/lincs-cell-painting/tree/master/metadata).

## Image-Based Profiling

We apply a unified, image-based profiling pipeline to all Drug Repurposing Hub Cell Painting plates.
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
