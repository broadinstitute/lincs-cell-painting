#!/usr/bin/env python
# coding: utf-8

# # Consensus signatures
# 
# Here, we generate consensus signatures for two batches of the LINCS Drug Repurposing Hub Cell Painting subset.
# Consensus signatures are level 5 data.
# 
# See the project [README.md](README.md) for more details.
# 
# This notebook generates eight files; one per plate normalization and consensus normalization strategy, with and without feature selection.
# 
# |Feature selection | Plate Normalization | Consensus Normalization | Consensus Suffix |
# |:---------------- | :------------------: | :------------------------: | -----------------: |
# | No  | DMSO | Median | `<BATCH>_consensus_median_dmso.csv.gz` |
# | No  | DMSO | MODZ | `<BATCH>_consensus_modz_dmso.csv.gz` |
# | No  | Whole Plate | Median | `<BATCH>_consensus_median.csv.gz` |
# | No  | Whole Plate | MODZ | `<BATCH>_consensus_modz.csv.gz` |
# | Yes | DMSO | Median | `<BATCH>_consensus_median_feature_select_dmso.csv.gz` |
# | Yes | DMSO | MODZ | `<BATCH>_consensus_modz_feature_select_dmso.csv.gz` |
# | Yes | Whole Plate | Median | `<BATCH>_consensus_median_feature_select.csv.gz` |
# | Yes | Whole Plate | MODZ | `<BATCH>_consensus_modz_feature_select.csv.gz` |
# 
# We generate eight files per batch.

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd

from pycytominer import aggregate, feature_select

from pycytominer import consensus
from pycytominer.cyto_utils import infer_cp_features, output


# In[3]:


# Set constants
batches = ["2016_04_01_a549_48hr_batch1", "2017_12_05_Batch2"]
operations = ["median", "modz"]
primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]

# Aggregating columns
replicate_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_cell_id",
    "Metadata_broad_sample",
    "Metadata_pert_well",
    "Metadata_mmoles_per_liter",
    "Metadata_dose_recode",
    "Metadata_time_point",
    "Metadata_moa",
    "Metadata_target",
]

# feature selection operations
feature_select_ops = [
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
]

# Output option
float_format = "%5g"
compression_options = {"method": "gzip", "mtime": 1}

# Two different blocklist feature sets (traditional and outlier)
commit = "838ac2eee8bee09b50a1ec8077201c01f7882c69"
traditional_blocklist_file = f"https://raw.githubusercontent.com/cytomining/pycytominer/{commit}/pycytominer/data/blocklist_features.txt"
outlier_blocklist_file = pathlib.Path("../utils/outlier_blocklist_features.txt")

full_blocklist_file = pathlib.Path("../utils/consensus_blocklist.txt")


# In[4]:


# Procees blocklist output
blocklist_df = pd.read_csv(traditional_blocklist_file)
outlier_blocklist_df = pd.read_csv(outlier_blocklist_file, names=["blocklist"])

full_blocklist_df = (
    pd.concat([blocklist_df, outlier_blocklist_df], axis="rows")
    .drop_duplicates()
    .reset_index(drop=True)
)

# Note, we need to output a file to be compatible with pycytominer
print(full_blocklist_df.shape)

full_blocklist_df.to_csv(full_blocklist_file, sep=",", index=False)


# In[5]:


# Set file information
file_bases = {
    "whole_plate": {
        "input_file_suffix": "_normalized.csv.gz",
        "output_file_suffix": ".csv.gz",
    },
    "dmso": {
        "input_file_suffix": "_normalized_dmso.csv.gz",
        "output_file_suffix": "_dmso.csv.gz",
    },
}

file_path_info = {"profile_dir": {}, "plate_dirs": {}, "plates": {}}
for batch in batches:
    file_path_info["profile_dir"][batch] = pathlib.Path("..", "profiles", batch)
    file_path_info["plate_dirs"][batch] = [
        x
        for x in file_path_info["profile_dir"][batch].iterdir()
        if x.name != ".DS_Store"
    ]
    file_path_info["plates"][batch] = [
        x.name for x in file_path_info["plate_dirs"][batch]
    ]
    print(len(file_path_info["plates"][batch]))

# The output directory is also the batch name
for batch in batches:
    pathlib.Path(batch).mkdir(exist_ok=True)


# ## Load and Process Data
# 
# We load data per plate, concatenate, add make minor modifications.
# 
# We perform this operation once per batch.

# In[6]:


# Load Data
all_profiles_dfs = {batch: {} for batch in batches}
cp_features = {batch: {} for batch in batches}

for batch in batches:
    for norm_strat, norm_file_base in file_bases.items():
        file_base = norm_file_base["input_file_suffix"]
        all_profiles_df = []
        for plate_dir in file_path_info["plate_dirs"][batch]:
            plate = plate_dir.name
            plate_file = plate_dir / f"{plate}{file_base}"
            plate_df = pd.read_csv(plate_file)
            all_profiles_df.append(plate_df)

        # Concatenate profiles
        all_profiles_df = pd.concat(all_profiles_df, axis="rows")

        # Add time metadata for batch 1 data
        if batch == "2016_04_01_a549_48hr_batch1":
            all_profiles_df = all_profiles_df.assign(Metadata_time_point="48H")

        # Recode missing MOA and target values to be "unknown"
        all_profiles_df.Metadata_moa = all_profiles_df.Metadata_moa.fillna("unknown")
        all_profiles_df.Metadata_target = all_profiles_df.Metadata_target.fillna(
            "unknown"
        )

        # Store concatenated data frame
        all_profiles_dfs[batch][norm_strat] = all_profiles_df

        # Determine every CellProfiler feature measured
        cp_features[batch][norm_strat] = infer_cp_features(
            all_profiles_dfs[batch][norm_strat]
        )

        # Clean up
        print(all_profiles_df.shape)
        del all_profiles_df


# ## Create consensus profiles
# 
# Create profiles with and without feature selection.
# 
# We generate two different consensus profiles for each of the normalization strategies, with and without feature selection.
# This generates eight different files _per batch_.

# In[7]:


all_consensus_dfs = {batch: {} for batch in batches}
for batch in batches:
    print(f"Now processing batch: {batch}")
    for norm_strat in file_bases:
        all_profiles_df = all_profiles_dfs[batch][norm_strat]
        cp_norm_features = cp_features[batch][norm_strat]

        consensus_profiles = {}
        for operation in operations:
            print(
                f"  Now calculating {operation} consensus for {norm_strat} normalization"
            )

            consensus_profiles[operation] = {}

            consensus_profiles[operation]["no_feat_select"] = consensus(
                profiles=all_profiles_df,
                replicate_columns=replicate_cols,
                operation=operation,
                features=cp_norm_features,
            )

            # How many DMSO profiles per well?
            print(
                f"  There are {consensus_profiles[operation]['no_feat_select'].shape[0]} {operation} consensus profiles for {norm_strat} normalization"
            )

            # Perform feature selection
            print(
                f"  Now feature selecting on {operation} consensus for {norm_strat} normalization"
            )
            consensus_profiles[operation]["feat_select"] = feature_select(
                profiles=consensus_profiles[operation]["no_feat_select"],
                features="infer",
                operation=feature_select_ops,
                blocklist_file=full_blocklist_file,
            )

            # How many features in feature selected profile?
            print(
                f"  There are {consensus_profiles[operation]['feat_select'].shape[1]} features in {operation} consensus profiles for {norm_strat} normalization"
            )

        all_consensus_dfs[batch][norm_strat] = consensus_profiles
    print("\n")


# ## Merge and output consensus signatures
# 
# Output with and without feature selection.

# In[8]:


for batch in batches:
    print(f"Now processing batch: {batch}")
    for norm_strat in file_bases:
        file_suffix = file_bases[norm_strat]["output_file_suffix"]
        for operation in operations:

            # No feature selection
            consensus_file = f"{batch}_consensus_{operation}{file_suffix}"
            consensus_file = pathlib.Path(batch, consensus_file)

            consensus_df = all_consensus_dfs[batch][norm_strat][operation][
                "no_feat_select"
            ]

            print(
                f"  Now Writing: Feature selection: No; Consensus Operation: {operation}; Norm Strategy: {norm_strat}"
            )
            print(f"  File: {consensus_file}")
            print(consensus_df.shape)

            output(
                df=consensus_df,
                output_filename=consensus_file,
                sep=",",
                float_format=float_format,
                compression_options=compression_options,
            )

            # With feature selection
            consensus_feat_df = all_consensus_dfs[batch][norm_strat][operation][
                "feat_select"
            ]

            consensus_feat_file = (
                f"{batch}_consensus_{operation}_feature_select{file_suffix}"
            )
            consensus_feat_file = pathlib.Path(batch, consensus_feat_file)

            print(
                f"  Now Writing: Feature selection: Yes; Consensus Operation: {operation}; Norm Strategy: {norm_strat}"
            )
            print(f"  File: {consensus_feat_file}")
            print(consensus_feat_df.shape)

            output(
                df=consensus_feat_df,
                output_filename=consensus_feat_file,
                sep=",",
                float_format=float_format,
                compression_options=compression_options,
            )
    print("\n")

