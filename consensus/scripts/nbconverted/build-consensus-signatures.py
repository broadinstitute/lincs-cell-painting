#!/usr/bin/env python
# coding: utf-8

# # Consensus Signatures
# 
# Here, we generate consensus signatures for the LINCS Drug Repurposing Hub Cell Painting subset.
# See the project [README.md](README.md) for more details.

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd

from pycytominer.aggregate import aggregate
from pycytominer.consensus import modz_base

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


# ## Load Profile Data Per Plate and Concatenate

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


# ## Create Consensus Profiles
# 
# We aggregate differently for DMSO and Compounds. Specifically, we aggregate by well for DMSO profiles in order to determine the extent of batch effect.

# In[7]:


# Aggregating columns
dmso_replicate_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_pert_well",
    "Metadata_mmoles_per_liter",
    "Metadata_dose_recode",
]

compound_replicate_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_broad_sample",
    "Metadata_mmoles_per_liter",
    "Metadata_dose_recode",
]


# In[8]:


all_consensus_dfs = {}
for norm_strat in file_bases:
    all_profiles_df = all_profiles_dfs[norm_strat]
    cp_norm_features = cp_features[norm_strat]

    # Isolate DMSO profiles
    dmso_profile_df = all_profiles_df.query(
        "Metadata_broad_sample == 'DMSO'"
    ).reset_index(drop=True)

    # How many DMSO profiles per well?
    print(
        f"There are {dmso_profile_df.shape[0]} DMSO profiles per well with {norm_strat} normalization"
    )

    # Form dmso consensus profiles with median and modz
    dmso_profiles = {}
    for operation in operations:
        print(
            f"Now calculating DMSO {operation} consensus for {norm_strat} normalization"
        )
        dmso_profiles[operation] = consensus_apply(
            dmso_profile_df,
            operation=operation,
            cp_features=cp_norm_features,
            replicate_cols=dmso_replicate_cols,
        )

    compound_profile_df = all_profiles_df.query(
        "Metadata_broad_sample != 'DMSO'"
    ).reset_index(drop=True)

    # Form compound consensus profiles with median and modz
    compound_profiles = {}
    for operation in operations:
        print(
            f"Now calculating Compound {operation} consensus for {norm_strat} normalization"
        )
        compound_profiles[operation] = consensus_apply(
            compound_profile_df,
            operation=operation,
            cp_features=cp_norm_features,
            replicate_cols=compound_replicate_cols,
        )

        compound_profiles[operation] = compound_profiles[operation].assign(
            Metadata_pert_well="collapsed"
        )

    consensus_profiles = {}
    for operation in operations:
        consensus_file = f"{batch}_consensus_{operation}.csv.gz"
        consensus_file = pathlib.Path(batch, consensus_file)

        dmso_df = dmso_profiles[operation]
        compound_df = compound_profiles[operation]

        consensus_profiles[operation] = (
            pd.concat([dmso_df, compound_df], axis="rows", sort=True)
            .reset_index(drop=True)
            .loc[:, dmso_replicate_cols + cp_norm_features]
        )

        # How many DMSO profiles per well?
        print(
            f"There are {consensus_profiles[operation].shape[0]} {operation} consensus profiles for {norm_strat} normalization"
        )

    all_consensus_dfs[norm_strat] = consensus_profiles


# ## Merge and Output Consensus Signatures

# In[9]:


for norm_strat in file_bases:
    file_suffix = file_bases[norm_strat]["output_file_suffix"]
    for operation in operations:
        consensus_file = f"{batch}_consensus_{operation}{file_suffix}"
        consensus_file = pathlib.Path(batch, consensus_file)

        consensus_df = all_consensus_dfs[norm_strat][operation]

        print(
            f"Now Writing: Consensus Operation: {operation}; Norm Strategy: {norm_strat}"
        )
        print(consensus_df.shape)

        consensus_df.to_csv(
            consensus_file, sep=",", compression="gzip", float_format="%5g", index=False
        )

