#!/bin/bash

python --version

# Reproduce the entire image-based profiling pipeline for batch 1 data
python profiling_pipeline.py --batch "2016_04_01_a549_48hr_batch1" --overwrite

# Apply the same profiling pipeline to process batch 2 data
python profiling_pipeline.py --batch "2017_12_05_Batch2" --plate_prefix "BR" --well_col "Metadata_Well" --plate_col "Metadata_Plate" --extract_cell_line --overwrite
