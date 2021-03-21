# Consensus Signatures

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

## Recoding Dose Information

The Drug Repurposing Hub collected data on 6 to 7 dose points per compound.
In general, most doses are very near the following 7 dose points (mmoles per liter):

> [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]

Therefore, to make it easier to filter by dose when comparing compounds, we first align the doses collected in the dataset to their nearest dose point above.
We then recode the dose points into ascending numerical levels and add a new metadata annotation `dose_recode` to the consensus signatures.

| Dose | Dose Recode |
| :--: | :---------: |
| 0 (DMSO) | 0 |
| ~0.04 | 1 |
| ~0.12 | 2 |
| ~0.37 | 3 |
| ~1.11 | 4 |
| ~3.33 | 5 |
| ~10 | 6 |
| ~20 | 7 |

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