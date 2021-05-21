# Consensus signatures

A consensus signature can be defined as a perturbation-specific summary profile acquired by aggregating replicate level information.
In the following notebook, we generate consensus profiles for each compound-dose pair.
We derive consensus signatures in two distinct ways:

1. Median Aggregation
2. Modified Z Score Aggregation (MODZ)

The first approach weights each replicate equally.
The second approach weights replicates by average similarity to other replicates.
See https://github.com/cytomining/pycytominer/issues/52 for more details and pointers to the MODZ calculation.

Also note that we form consensus profiles for two different normalization strategies: Whole plate z scoring and DMSO z scoring.
Therefore, we generate a total of four consensus signature files.

## DMSO vs. Compound Profiles

Note we generated per-well DMSO consensus signatures and per compound-dose pair consensus signatures for compounds.
The per-well DMSO profiles can help to assess plate-associated batch effects.

## Reproduce Pipeline

The pipeline can be reproduced by executing the following:

```bash
# Make sure conda environment is activated
conda activate lincs

# Reproduce the pipeline for producing bulk signatures
ipython scripts/nbconverted/build-consensus-signatures.py
```

`scripts/nbconverted/*.py` were created from the Jupyter notebooks in this folder, like this:

```sh
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
```
