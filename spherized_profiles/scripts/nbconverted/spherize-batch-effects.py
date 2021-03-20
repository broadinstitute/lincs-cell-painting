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
    "blacklist",
    "drop_outliers"
]

na_cut = 0
corr_threshold = 0.95
outlier_cutoff = 60

output_dir = "profiles"


# In[3]:


for batch in batches:
    for suffix in suffixes:
        output_file = pathlib.Path(
            f"{output_dir}/{batch}_dmso_spherized_profiles_with_input_normalized_by_{suffix}.csv.gz"
        )
        print(f"Now processing {output_file}...")

        profile_df = pd.concat([pd.read_csv(x) for x in files[batch][suffix]]).reset_index(drop=True)

        # Perform feature selection
        profile_df = feature_select(
            profiles=profile_df,
            operation=feature_select_ops,
            na_cutoff=0,
            corr_threshold=corr_threshold,
            outlier_cutoff=outlier_cutoff
        )

        print(profile_df.shape)
        profile_df.head()

        spherize_df = normalize(
            profiles=profile_df,
            features="infer",
            meta_features="infer",
            samples="Metadata_broad_sample == 'DMSO'",
            method="whiten",
        )

        print(spherize_df.shape)
        spherize_df.head()

        output(
            df=spherize_df,
            output_filename=output_file
        )

