# Image-Based profiling

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

## Data Levels

We include two batches of Cell Painting data in this repository: `2016_04_01_a549_48hr_batch1` and `2017_12_05_Batch2`.

### CellProfilier-derived profiles

For each batch, we include:

| Data level | Description | File format | Included in this repo? |
| :--------- | :---------- | :---------- | :-------------------- |
| Level 1 | Cell images | `.tif` | No^ |
| Level 2 | Single cell profiles | `.sqlite` | No^ |
| Level 3 | Aggregated profiles with metadata | `.csv.gz` | Yes |
| Level 4a | Normalized profiles with metadata | `.csv.gz` | Yes |
| Level 4b | Normalized and feature selected profiles with metadata | `.csv.gz` | Yes |
| Level 5 | Consensus perturbation profiles | `.csv.gz` | Yes |

Importantly, we include files for _two_ different types of normalization: Whole-plate normalization, and DMSO-specific normalization.
See [`profile_cells.py`](profile_cells.py) for more details.

#### Batch corrected profiles

We use a spherize (a.k.a. whiten) transform to adjust for plate position effects.
The spherize transform adjusts for plate position effects by transforming the profile data such that the DMSO profiles are left with an identity covariance matrix.
See [`spherize-batch-effects.ipynb`](spherized_profiles/spherize-batch-effects.ipynb) for implementation details.

For each batch we include four different spherized profiles.
These data include all level 4b profiles for every batch.

| Batch | Input data | Spherized output file |
| :---: | :--------: | :-------------------: |
| 2016_04_01_a549_48hr_batch1 | DMSO normalized | 2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.csv.gz |
| 2016_04_01_a549_48hr_batch1 | Whole plate normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz |
| 2017_12_05_Batch2 | DMSO normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.csv.gz |
| 2017_12_05_Batch2 | Whole plate normalized | 2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz |

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

## Critical details

There are several critical details that are important for understanding data generation and processing.
See [`profile_cells.py`](profile_cells.py) for more details about the specific processing steps and decisions.
