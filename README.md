# Processed Cell Painting Data for the LINCS Drug Repurposing Project

The repository stores data and data processing scripts for **a subset** of the LINCS Drug Repurposing Project.
The specific subset include all Cell Painting profiles for ~1,500 compounds.

For a specific list of compounds tested, see [`metadata`](https://github.com/broadinstitute/lincs-cell-painting/tree/master/metadata).

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
