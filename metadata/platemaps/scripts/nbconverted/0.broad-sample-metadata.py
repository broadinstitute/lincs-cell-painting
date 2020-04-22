#!/usr/bin/env python
# coding: utf-8

# ## Extract all Metadata for Broad IDs Assayed with Cell Painting

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pandas as pd


# In[3]:


all_platemaps = list()

platemap_dir = "platemap"
for platemap_file in os.listdir(platemap_dir):

    # Load platemap
    platemap = platemap_file.strip(".txt")
    platemap_df = pd.read_csv(os.path.join(platemap_dir, platemap_file), sep="\t")

    assert platemap == platemap_df.plate_map_name.unique()[0]

    # Process platemap
    platemap_df = platemap_df.assign(
        broad_id=platemap_df.broad_sample.str.extract(r"(BRD[-N][A-Z0-9]+)")
    )

    platemap_df = (
        platemap_df.loc[:, ["broad_sample", "broad_id", "plate_map_name", "solvent"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    all_platemaps.append(platemap_df)


# In[4]:


# Combine platemap info
all_platemap_df = pd.concat(all_platemaps, axis="rows").drop_duplicates()

# Output file
output_file = "broad_sample_info.tsv"
all_platemap_df.to_csv(output_file, index=False, sep="\t")

print(all_platemap_df.shape)
all_platemap_df.head()

