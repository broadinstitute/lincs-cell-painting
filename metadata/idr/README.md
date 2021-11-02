# IDR submission

_TODO: Batch 2 files are currently missing_

This folder contains the following files

- input files, R notebooks (one per batch of data) for generating the library metadata file for submitting the imaging data to IDR, and the output of the notebook
- manually generated study files `idr0000-study_batch_1.csv` and `idr0000-study_batch_2.csv`

The input files for the notebook are as follows:

- In this folder:
  - `lincs_plates_final_batch_1.txt` and `lincs_plates_final_batch_2.txt` - list of full plate names
  - `idr0000-screenA-library.txt` - template library file taken from https://github.com/IDR/idr0000-lastname-example/archive/master.zip
- At other locations in this repository: 
  - `repurposing_drugs_20200324.txt` - repurposing drug annotations
  - `repurposing_samples_20200324.txt` - repurposing sample annotations
  - `2016_04_01_a549_48hr_batch1_metadata_cell_count_summary.tsv.gz` and `2017_12_05_Batch2_metadata_cell_count_summary.tsv.gz` - cell counts, but we use it here as an easy-to-read source of raw screen metadata

The output file  is:

- `idr0000-screenA-library_batch_1.txt.gz` and `idr0000-screenA-library_batch_2.txt.gz` – using `idr0000-screenA-library.txt` as the template

