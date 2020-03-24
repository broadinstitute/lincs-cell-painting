#!/usr/bin/env python
# coding: utf-8

# # Merging Compound files distributed by CLUE
# 
# We consolidate drug and sample resources information into a single file for easier downstream processing.
# 
# The data were originally retrieved from https://clue.io/repurposing.
# See [`clue/README.md`](clue/README.md) for more details.

# In[1]:


import os
import numpy as np
import pandas as pd


# ## Load Data

# In[2]:


data_dir = "clue"
date = "20200322"


# In[3]:


drug_file = os.path.join(data_dir, "repurposing_drugs_{}.txt".format(date))
drug_df = pd.read_csv(drug_file, encoding = "ISO-8859-1", sep='\t', skiprows=9)

print(drug_df.shape)
drug_df.head(2)


# In[4]:


sample_file = os.path.join(data_dir, "repurposing_samples_{}.txt".format(date))
sample_df = pd.read_csv(sample_file, encoding = "ISO-8859-1", sep='\t', skiprows=9)

print(sample_df.shape)
sample_df.head(2)


# ## Checking for `pert_iname` Discrepancies

# In[5]:


# Assert that all pert_inames exist in both resources
assert len(set(drug_df.pert_iname.values).difference(set(sample_df.pert_iname))) == 0
assert len(set(sample_df.pert_iname.values).difference(set(drug_df.pert_iname))) == 0


# ## Merge the Samples and Drugs data

# In[6]:


combined_df = (
    drug_df.merge(
        sample_df,
        on="pert_iname",
        how="inner"
    )
    .reset_index(drop=True)
)

# Move broad_id to first column
col_order = combined_df.columns.tolist()
col_order.insert(0, col_order.pop(col_order.index('broad_id')))
combined_df = combined_df.loc[:, col_order]

# Output to file
output_file = "repurposing_info"
combined_df.to_csv("{}.tsv".format(output_file), sep='\t', index=False)

print(combined_df.shape)
combined_df.head()


# ## Create a "Long" version where we split MOA and Target delimiters
# 
# Certain compounds have multiple MOA classes and targets that are delimited by pipes (`|`).
# Each MOA class and target can be considered to have equal support (see https://github.com/broadinstitute/lincs-cell-painting/issues/5).
# 
# Split the combined data on both MOA and target along each pipe and elongate the table.
# This is done to reduce computational burden of multiple downstream analyses performing the same splits.

# In[7]:


# The splitting strategy does not work with missing values
# Add a dummy variable, that will be replaced downstream
combined_df.moa = combined_df.moa.fillna("replace_with_na")
combined_df.target = combined_df.target.fillna("replace_with_na")


# In[8]:


# Make sure the original index is preserved
split_col_index = "{}_index".format(output_file)


# In[9]:


moa_split_df = (
    pd.DataFrame(combined_df.moa.str.split("|").tolist(), index=combined_df.index)
    .stack()
    .reset_index()
)
moa_split_df.columns = [split_col_index, "_", "moa_unique"]

print(moa_split_df.shape)
moa_split_df.head()


# In[10]:


target_split_df = (
    pd.DataFrame(combined_df.target.str.split("|").tolist(), index=combined_df.index)
    .stack()
    .reset_index()
)

target_split_df.columns = [split_col_index, "_", "target_unique"]

print(target_split_df.shape)
target_split_df.head()


# In[11]:


long_combined_df = (
    combined_df
    .merge(
        moa_split_df.loc[:, [split_col_index, "moa_unique"]],
        left_index=True,
        right_on=split_col_index,
        how="left"
    )
    .merge(
        target_split_df.loc[:, [split_col_index, "target_unique"]],
        on=split_col_index,
        how="left"
    )
    .reset_index(drop=True)
)

# Put back missing values
long_combined_df.loc[long_combined_df.moa == "replace_with_na", "moa"] = np.nan
long_combined_df.loc[long_combined_df.moa_unique == "replace_with_na", "moa_unique"] = np.nan
long_combined_df.loc[long_combined_df.target == "replace_with_na", "target"] = np.nan
long_combined_df.loc[long_combined_df.target_unique == "replace_with_na", "target_unique"] = np.nan

# Output to file
output_file = "repurposing_info_long.tsv"
long_combined_df.to_csv(output_file, sep='\t', index=False)

print(long_combined_df.shape)
long_combined_df.head()

