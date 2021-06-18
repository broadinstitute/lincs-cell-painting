# Image-based profiling

Image-based profiling represents a series of data processing steps that turn image-based readouts into more manageable data matrices for downstream analyses ([Caicedo et al. 2017](https://doi.org/10.1038/nmeth.4397)).
Typically, you derive image-based readouts using software, like CellProfiler ([McQuin et al. 2018](https://doi.org/10.1371/journal.pbio.2005970)), that segment cells and extract so-called hand-engineered single cell morphology measurements.
In this folder, we process the CellProfiler derived morphology features for the LINCS Cell Painting dataset using [pycytominer](https://github.com/cytomining/pycytominer) - a tool enabling reproducible image-based profiling.

Specifically, we include:

1. Data processing scripts to perform the full unified, image-based profiling pipeline
2. Processed data for each Cell Painting plate (for several "data levels")
3. Instructions on how to reproduce the profiling pipeline

## Workflow

![Cytominer workflow](media/cytominer_workflow.png)

Note here that we do not include the intermediate step of generating `.sqlite` files per plate using a tool called [cytominer-database](https://github.com/cytomining/cytominer-database).
This repository and workflow begins after we applied cytominer-database.

## Data levels

We include two batches of Cell Painting data in this repository: `2016_04_01_a549_48hr_batch1` and `2017_12_05_Batch2`.

### CellProfilier-derived profiles

For each batch, we include:

| Data level | Description | File format | Included in this repo? | Version control |
| :--------- | :---------- | :---------- | :--------------------- | :-------------- |
| Level 1 | Cell images | `.tif` | No^ | NA |
| Level 2 | Single cell profiles | `.sqlite` | No^ | NA |
| Level 3 | Aggregated profiles with metadata | `.csv.gz` | Yes | dvc |
| Level 4a | Normalized profiles with metadata | `.csv.gz` | Yes | dvc |
| Level 4b | Normalized and feature selected profiles with metadata | `.csv.gz` | Yes | dvc |
| Level 5 | Consensus perturbation profiles | `.csv.gz` | Yes | git lfs |

Importantly, we include files for _two_ different types of normalization: Whole-plate normalization, and DMSO-specific normalization.
See [`profile_cells.py`](profile_cells.py) for more details.

> Note: If you use normalized profiles without feature selection, you may need to remove additional outlier features.
See https://github.com/broadinstitute/lincs-cell-painting/issues/65 for more details and a list of features.

#### Batch corrected profiles

We use a spherize (a.k.a. whiten) transform to adjust for plate position effects.
The spherize transform adjusts for plate position effects by transforming the profile data such that the DMSO profiles are left with an identity covariance matrix.
See [`spherize-batch-effects.ipynb`](spherized_profiles/spherize-batch-effects.ipynb) for implementation details.

For each batch we include four different spherized profiles.
These data include all level 4b profiles for every batch.

| Batch | Input data | Spherized output file | Version control |
| :---: | :--------: | :-------------------: | :-------------- |
| 2016_04_01_a549_48hr_batch1 | DMSO normalized | 2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.csv.gz | git lfs |
| 2016_04_01_a549_48hr_batch1 | Whole plate normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz | git lfs |
| 2017_12_05_Batch2 | DMSO normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.csv.gz | git lfs |
| 2017_12_05_Batch2 | Whole plate normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz | git lfs |

### Data access

We use a combination of [`git lfs`](https://git-lfs.github.com/) and [`dvc`](https://dvc.org/) to version all the data in this repository (see tables above for a breakdown).
In order to access the data locally, you must first install both services. 

```bash
# Download consensus profiles and spherized data
git lfs pull

# Download individual plate data
dvc pull
```

### DeepProfiler-derived profiles

TBD

## Reproduce pipeline

After activating the conda environment with:

```bash
conda activate lincs
```

you can reproduce the pipeline by simply executing the following:

```bash
# Make sure you are in the profiles/ directory
./run.sh
```

## Recoding dose information

The Drug Repurposing Hub collected data on 6 to 7 dose points per compound.
In general, most doses are very near the following 7 dose points (mmoles per liter):

> [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]

Therefore, to make it easier to filter by dose when comparing compounds, we first align the doses collected in the dataset to their nearest dose point above.
We then recode the dose points into ascending numerical levels and add a new metadata annotation `Metadata_dose_recode` to all profiles.

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

## Critical details

There are several critical details that are important for understanding data generation and processing.
See [`profile_cells.py`](profile_cells.py) for more details about the specific processing steps and decisions.
