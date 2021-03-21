#!/usr/bin/env python
# coding: utf-8

# # Consensus Signatures
# 
# Here, we generate consensus signatures for the LINCS Drug Repurposing Hub Cell Painting subset.
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

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd

from pycytominer.aggregate import aggregate
from pycytominer.consensus import modz_base
from pycytominer.feature_select import feature_select
from pycytominer.cyto_utils import infer_cp_features


# In[3]:


def recode_dose(x, doses, return_level=False):
    closest_index = np.argmin([np.abs(dose - x) for dose in doses])
    if np.isnan(x):
        return 0
    if return_level:
        return closest_index + 1
    else:
        return doses[closest_index]


def consensus_apply(df, operation, cp_features, replicate_cols):
    if operation == "modz":
        consensus_df = (
            df.groupby(replicate_cols)
            .apply(lambda x: modz_base(x.loc[:, cp_features]))
            .reset_index()
        )
    elif operation == "median":
        consensus_df = aggregate(
            df, operation="median", features="infer", strata=replicate_cols
        )
    return consensus_df


# In[4]:


# Set constants
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
operations = ["median", "modz"]
batch = "2016_04_01_a549_48hr_batch1"
primary_dose_mapping = [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]

# Set file paths
profile_dir = pathlib.Path("..", "profiles", batch)
plate_dirs = [x for x in profile_dir.iterdir() if x.name != ".DS_Store"]
plates = [x.name for x in plate_dirs]
print(len(plates))


# In[5]:


# The output directory is also the batch name
pathlib.Path(batch).mkdir(exist_ok=True)


# ## Load and Process Data
# 
# We load data per plate, concatenate, and recode dose information

# In[6]:


# Load Data
all_profiles_dfs = {}
cp_features = {}
for norm_strat, norm_file_base in file_bases.items():
    file_base = norm_file_base["input_file_suffix"]
    all_profiles_df = []
    for plate_dir in plate_dirs:
        plate = plate_dir.name
        plate_file = plate_dir / f"{plate}{file_base}"
        plate_df = pd.read_csv(plate_file)
        all_profiles_df.append(plate_df)

    # Concatenate profiles
    all_profiles_df = pd.concat(all_profiles_df, axis="rows")

    # Recode dose
    all_profiles_df = all_profiles_df.assign(
        Metadata_dose_recode=(
            all_profiles_df.Metadata_mmoles_per_liter.apply(
                lambda x: recode_dose(x, primary_dose_mapping, return_level=True)
            )
        )
    )

    # Make sure DMSO profiles recieve a zero dose level
    all_profiles_df.loc[
        all_profiles_df.Metadata_broad_sample == "DMSO", "Metadata_dose_recode"
    ] = 0

    # Store concatenated data frame
    all_profiles_dfs[norm_strat] = all_profiles_df

    # Determine every CellProfiler feature measured
    cp_features[norm_strat] = infer_cp_features(all_profiles_dfs[norm_strat])

    # Clean up
    print(all_profiles_df.shape)
    del all_profiles_df


# ## Create Consensus Profiles, with and without feature selection
# 
# We generate two different consensus profiles for each of the normalization strategies, with and without feature selection. This generates eight different files.

# In[7]:


# Aggregating columns
replicate_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_pert_well",
    "Metadata_mmoles_per_liter",
    "Metadata_dose_recode",
    "Metadata_moa",
    "Metadata_target",
]


# In[8]:


# feature selection operations
feature_select_ops = [
    "drop_na_columns",
    "variance_threshold",
    "correlation_threshold",
    "blacklist",
]

all_consensus_dfs = {}
for norm_strat in file_bases:
    all_profiles_df = all_profiles_dfs[norm_strat]
    cp_norm_features = cp_features[norm_strat]

    consensus_profiles = {}
    for operation in operations:
        print(f"Now calculating {operation} consensus for {norm_strat} normalization")

        consensus_profiles[operation] = {}

        consensus_profiles[operation]["no_feat_select"] = consensus_apply(
            all_profiles_df,
            operation=operation,
            cp_features=cp_norm_features,
            replicate_cols=replicate_cols,
        )

        # How many DMSO profiles per well?
        print(
            f"There are {consensus_profiles[operation]['no_feat_select'].shape[0]} {operation} consensus profiles for {norm_strat} normalization"
        )

        # feature selection
        print(
            f"Now feature selecting on {operation} consensus for {norm_strat} normalization"
        )

        consensus_profiles[operation]["feat_select"] = feature_select(
            profiles=consensus_profiles[operation]["no_feat_select"],
            features="infer",
            operation=feature_select_ops,
        )

        # How many features in feature selected profile?
        print(
            f"There are {consensus_profiles[operation]['feat_select'].shape[1]} features in {operation} consensus profiles for {norm_strat} normalization"
        )

    all_consensus_dfs[norm_strat] = consensus_profiles


# ## Merge and Output Consensus Signatures, with and without feature selection

# In[9]:


float_format = "%5g"
compression = "gzip"

for norm_strat in file_bases:
    file_suffix = file_bases[norm_strat]["output_file_suffix"]
    for operation in operations:

        # No feature selection
        consensus_file = f"{batch}_consensus_{operation}{file_suffix}"
        consensus_file = pathlib.Path(batch, consensus_file)

        consensus_df = all_consensus_dfs[norm_strat][operation]["no_feat_select"]

        print(
            f"Now Writing: Feature selection: No; Consensus Operation: {operation}; Norm Strategy: {norm_strat}\nFile: {consensus_file}"
        )
        print(consensus_df.shape)

        consensus_df.to_csv(
            consensus_file,
            sep=",",
            compression=compression,
            float_format=float_format,
            index=False,
        )

        # With feature selection
        consensus_feat_df = all_consensus_dfs[norm_strat][operation]["feat_select"]

        consensus_feat_file = (
            f"{batch}_consensus_{operation}_feature_select{file_suffix}"
        )
        consensus_feat_file = pathlib.Path(batch, consensus_feat_file)

        print(
            f"Now Writing: Feature selection: Yes; Consensus Operation: {operation}; Norm Strategy: {norm_strat}\nFile: {consensus_feat_file}"
        )
        print(consensus_feat_df.shape)

        consensus_feat_df.to_csv(
            consensus_feat_file,
            sep=",",
            compression=compression,
            float_format=float_format,
            index=False,
        )

