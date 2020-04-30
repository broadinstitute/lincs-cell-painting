# Image-Based Profiling

Image-based profiling represents a series of data processing steps that turn Image-based readouts into more manageable data matrices for downstream analyses ([Caicedo et al. 2017](https://doi.org/10.1038/nmeth.4397)).
Typically, the image-based readouts are derived from CellProfiler ([McQuin et al. 2018](https://doi.org/10.1371/journal.pbio.2005970)) and represent single cell morphology measurements.
In this folder, we process the CellProfiler derived morphology features using [pycytominer](https://github.com/cytomining/pycytominer) - a tool enabling reproducible image-based profiling.
Specifically, we include:

1. Data processing scripts to perform the full unified, image-based profiling pipeline
2. Processed data for each Cell Painting plate (for several "data levels")
3. Instructions on how to reproduce the profiling pipeline

## Workflow

![Cytominer Workflow](media/cytominer_workflow.png)

Note here that we do not include the intermediate step of generating `.sqlite` files per plate using a tool called [cytominer-database](https://github.com/cytomining/cytominer-database).
This repository and workflow begins after we applied cytominer-database.

## Data Levels

### CellProfilier-derived Profiles

| Data Level | Description | File Format | Included in this Repo |
| :--------- | :---------- | :---------- | :-------------------- |
| Level 1 | Cell Images | `.tif` | No^ |
| Level 2 | Single Cell Profiles | `.sqlite` | No^ |
| Level 3 | Aggregated Profiles with Metadata | `.csv.gz` | Yes |
| Level 4a | Normalized Profiles with Metadata | `.csv.gz` | Yes |
| Level 4b | Normalized and Feature Selected Profiles with Metadata | `.csv.gz` | Yes |
| Level 5 | Consensus Perturbation Profiles | `.csv.gz` | Yes |

^ Note that these files are being prepared

### DeepProfiler-derived Profiles

TBD

## Reproduce Pipeline

The pipeline can be reproduced by simply executing the following:

```bash
# Make sure conda environment is activated
conda activate lincs

# Reproduce the entire image-based profiling pipeline for CellProfiler derived features
python profiling_pipeline.py
```

## Critical Details

There are several critical details that are important for understanding data generation and processing.
See [`profile.py`](profile.py) for more details about the specific processing steps and decisions.
