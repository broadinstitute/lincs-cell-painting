#!/usr/bin/env python
# coding: utf-8

# ## Adjust batch effects with a spherize transform
# 
# Here, we load in all normalized profiles (level 4a) data across all plates and apply a spherize transform using the DMSO profiles as the background distribution.
# 
# We've previously observed that sphering (aka whitening) the data successfully adjusts for technical artifacts induced by batch to batch variation and plate position effects.

# In[1]:


import os
import pathlib
import subprocess
import pandas as pd

from pycytominer import normalize, feature_select
from pycytominer.cyto_utils import output, infer_cp_features


# In[2]:


input_dir = pathlib.Path("../profiles/")
batches = ["2016_04_01_a549_48hr_batch1", "2017_12_05_Batch2"]

suffixes = {
    "whole_plate": "_normalized.csv.gz",
    "dmso": "_normalized_dmso.csv.gz"
}

plates = {
    batch: [x.name for x in pathlib.Path(f"{input_dir}/{batch}").iterdir() if ".DS_Store" not in x.name]
    for batch in batches
}

files = {
    batch: {
        suffix: [pathlib.Path(f"{input_dir}/{batch}/{x}/{x}{suffixes[suffix]}") for x in plates[batch]]
        for suffix in suffixes
    }
    for batch in batches
}

feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "drop_na_columns",
    "blocklist",
]

na_cut = 0
corr_threshold = 0.95
output_dir = "profiles"

full_blocklist_file = pathlib.Path("../utils/consensus_blocklist.txt")


# In[3]:


for batch in batches:
    for suffix in suffixes:
        output_file = pathlib.Path(
            f"{output_dir}/{batch}_dmso_spherized_profiles_with_input_normalized_by_{suffix}.csv.gz"
        )
        print(f"Now processing {output_file}...")

        profile_df = pd.concat([pd.read_csv(x) for x in files[batch][suffix]]).reset_index(drop=True)
        print(profile_df.shape)
        
        # Step 1: Perform feature selection
        if batch == "2017_12_05_Batch2":
            profile_df = (
                profile_df
                .groupby(["Metadata_cell_line", "Metadata_time_point"])
                .apply(
                    lambda x: feature_select(
                        profiles=x,
                        operation=feature_select_ops,
                        na_cutoff=na_cut,
                        corr_threshold=corr_threshold,
                        blocklist_file=full_blocklist_file
                    )
                )
            )
            
            # Drop features that weren't selected in the grouped splits
            profile_df = feature_select(
                profiles=profile_df,
                operation="drop_na_columns",
                na_cutoff=na_cut
            )
        else:
            profile_df = feature_select(
                profiles=profile_df,
                operation=feature_select_ops,
                na_cutoff=na_cut,
                corr_threshold=corr_threshold,
                blocklist_file=full_blocklist_file
            )

        # Step 2: Spherize transform
        if batch == "2017_12_05_Batch2":
            spherize_df = (
                profile_df
                .groupby(["Metadata_cell_line", "Metadata_time_point"])
                .apply(
                    lambda x: normalize(
                        profiles=x,
                        features="infer",
                        meta_features="infer",
                        samples="Metadata_broad_sample == 'DMSO'",
                        method="spherize"
                    )
                )
            )
        else:
            spherize_df = normalize(
                profiles=profile_df,
                features="infer",
                meta_features="infer",
                samples="Metadata_broad_sample == 'DMSO'",
                method="spherize"
            )

        print(spherize_df.shape)
        spherize_df.head()

        # Step 3: Output profiles
        output(
            df=spherize_df,
            output_filename=output_file
        )

