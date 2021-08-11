#!/usr/bin/env python
# coding: utf-8

# ## Generate consensus signatures of spherized profiles
# 
# In `spherize-batch-effects.ipynb`, we performed feature selection and spherized level 4 profiles.
# Here, we acquire consensus signatures for these spherized data.

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd

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

features = "infer"

# Output options
float_format = "%5g"
compression_options = {"method": "gzip", "mtime": 1}
output_dir = pathlib.Path("consensus")

# Two different blocklist feature sets (traditional and outlier)
commit = "838ac2eee8bee09b50a1ec8077201c01f7882c69"
traditional_blocklist_file = f"https://raw.githubusercontent.com/cytomining/pycytominer/{commit}/pycytominer/data/blocklist_features.txt"
outlier_blocklist_file = pathlib.Path("../utils/outlier_blocklist_features.txt")

full_blocklist_file = pathlib.Path("../utils/consensus_blocklist.txt")


# In[4]:


# Establish input files
profile_dir = pathlib.Path("profiles")
spherized_string = "_dmso_spherized_profiles_with_input_normalized_by_"
spherized_suffixes = ["dmso", "whole_plate"]

profile_files = {
    batch: {
        suffix: pathlib.Path(f"{profile_dir}/{batch}{spherized_string}{suffix}.csv.gz")
        for suffix in spherized_suffixes
    }
    for batch in batches
}

profile_files


# In[5]:


all_consensus_dfs = {batch: {} for batch in batches}
for batch in batches:
    print(f"Now processing batch: {batch}")
    batch_files = profile_files[batch]
    for norm_strat in batch_files:
        spherized_file = batch_files[norm_strat]
        print(f"  Now forming consensus signature for: {spherized_file}")

        spherized_df = pd.read_csv(spherized_file, low_memory=False)
        print(spherized_df.shape)

        # Recode missing MOA and target values to be "unknown"
        spherized_df.Metadata_moa = spherized_df.Metadata_moa.fillna("unknown")
        spherized_df.Metadata_target = spherized_df.Metadata_target.fillna("unknown")

        # Set a timepoint variable only for batch 1
        if batch == "2016_04_01_a549_48hr_batch1":
            spherized_df = spherized_df.assign(Metadata_time_point="48H")

        for operation in operations:
            output_file = pathlib.Path(
                f"{output_dir}/{batch}{spherized_string}{norm_strat}_consensus_{operation}.csv.gz"
            )
            print(f"    with consensus operation: {operation}")

            spherized_consensus_df = consensus(
                profiles=spherized_df,
                replicate_columns=replicate_cols,
                operation=operation,
                features=features,
            )
            print(spherized_consensus_df.shape)

            output(
                df=spherized_consensus_df,
                output_filename=output_file,
                sep=",",
                float_format=float_format,
                compression_options=compression_options,
            )
            print("    Done.")

    print("Batch done.\n")

