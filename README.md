# LINCS Cell Painting profile data repository

The Library of Integrated Network-Based Cellular Signatures (LINCS) Project aims to create publicly available resources to characterize how cells respond to perturbation.
This repository stores Cell Painting readouts and associated data-processing pipelines for the LINCS Cell Painting dataset.

In this project, the [Connectivity Map](https://clue.io/team) team perturbed A549 cells with 1,571 compounds across 6 doses in 5 technical replicates.
The data represent **a subset** of the [Broad Drug Repurposing Hub](https://clue.io/repurposing#home) collection of compounds.

We refer to this dataset as `LINCS Pilot 1`.
We also include data for the second batch of LINCS Cell Painting data, which we refer to as `LKCP`.

For a specific list of compounds tested, see [`metadata`](https://github.com/broadinstitute/lincs-cell-painting/tree/master/metadata).
You can interactively explore information about the compounds in the [CLUE Repurposing app](https://clue.io/repurposing-app).

The [Morphology Connectivity Hub](https://clue.io/morphology) is the primary source of this dataset.

## Image-based profiling

We apply a unified, image-based profiling pipeline to all 136 384-well plates from `LINCS Pilot 1`, and all 135 384-well plates from `LKCP`.
We use [pycytominer](https://github.com/cytomining/pycytominer) as the primary tool for image-based profiling.

We process and store level 3 to level 5 profiles in the [profiles/](profiles/) directory.
Furthermore, spherized and conensus profiles can be found in their relevant folders.

See [`profiles/README.md`](profiles/README.md) for more details and for instructions on how to reproduce the pipeline.
For further details about image-based profiling in general, please refer to [Caicedo et al. 2017](https://doi.org/10.1038/nmeth.4397).

## Computational environment

We use [conda](https://docs.conda.io/en/latest/) to manage the computational environment.

To install conda see [instructions](https://docs.conda.io/en/latest/miniconda.html).

We recommend installing conda by downloading and executing the `.sh` file and accepting defaults.

After installing conda, execute the following to install and navigate to the environment:

```bash
# First, install the `lincs` conda environment
conda env create --force --file environment.yml

# If you had already installed this environment and now want to update it
conda env update --file environment.yml --prune

# Then, activate the environment and you're all set!
conda activate lincs
```

Also note that when contributing to the repository, make sure to add any new package in the `environment.yml` file.

## License

We use a dual license in this repository.
We license the source code as [BSD 3-Clause](LICENSE_BSD3.md), and license the data, results, and figures as [CC0 1.0](LICENSE_CC0.md).

## Citation

If you use these data or software, please cite our [Zenodo archive](https://zenodo.org/record/5008187):

> Natoli, Ted, Way, Gregory, Lu, Xiaodong, Logan, David, Alimova, Maria, Hartland, Kate, Golub, Todd, Carpenter, Anne, Singh, Shantanu, Subramanian, Aravind. (2021). broadinstitute/lincs-cell-painting: Full release of LINCS Cell Painting dataset (Version v1) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.5008187
