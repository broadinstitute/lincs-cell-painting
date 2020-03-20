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
date = "20180907"


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


assert len(set(drug_df.pert_iname.values).difference(set(sample_df.pert_iname))) == 0
set(sample_df.pert_iname.values).difference(set(drug_df.pert_iname))


# Two perturbation names (`pert_iname`) are inconsistent.
# Work towards reconciliation.

# ### `YM-298198-desmethyl` 

# In[6]:


sample_df.loc[sample_df.pert_iname.str.contains("298198"), :]


# In[7]:


ym_drug = drug_df.loc[drug_df.pert_iname.str.contains("298198"), :].reset_index(drop=True)
ym_drug.loc[0, "pert_iname"] = "YM-298198-desmethyl"
ym_drug


# In[8]:


# solve the YM-298198 problem
drug_df = pd.concat([drug_df, ym_drug], axis="rows").reset_index(drop=True)


# #### `YM-298198-desmethyl` is absent in the drug data
# 
# [`YM-298198-desmethyl`](https://www.tocris.com/products/desmethyl-ym-298198_2447) is a derivative of [`YM-298198`](https://www.tocris.com/products/ym-298198-hydrochloride_2448), and therefore has a different structure. However, their MOA and target is the same.

# ### `golgicide-A`

# In[9]:


sample_df.loc[sample_df.pert_iname.str.contains("golgicide"), :]


# In[10]:


drug_df.loc[drug_df.pert_iname.str.contains("golgicide"), :]


# #### `golgicide` only differs by a capitalization and is equivalent

# In[11]:


# Solve the golgicide problem
sample_df.loc[sample_df.pert_iname.str.contains("golgicide"), "pert_iname"] = "golgicide-a"


# In[12]:


# Now, assert that there are no differences
assert len(set(sample_df.pert_iname.values).difference(set(drug_df.pert_iname))) == 0


# ## Merge the Samples and Drugs data

# In[13]:


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
combined_df.to_csv(".tsv".format(output_file), sep='\t', index=False)

print(combined_df.shape)
combined_df.head()


# ## Create a "Long" version where we split MOA and Target delimiters
# 
# Certain compounds have multiple MOA classes and targets that are delimited by pipes (`|`).
# Each MOA class and target can be considered to have equal support (see https://github.com/broadinstitute/lincs-cell-painting/issues/5).
# 
# Split the combined data on both MOA and target along each pipe and elongate the table.
# This is done to reduce computational burden of multiple downstream analyses performing the same splits.

# In[14]:


# The splitting strategy does not work with missing values
# Add a dummy variable, that will be replaced downstream
combined_df.moa = combined_df.moa.fillna("replace_with_na")
combined_df.target = combined_df.target.fillna("replace_with_na")


# In[15]:


# Make sure the original index is preserved
split_col_index = "{}_index".format(output_file)


# In[16]:


moa_split_df = (
    pd.DataFrame(combined_df.moa.str.split("|").tolist(), index=combined_df.index)
    .stack()
    .reset_index()
)
moa_split_df.columns = [split_col_index, "_", "moa_unique"]

print(moa_split_df.shape)
moa_split_df.head()


# In[17]:


target_split_df = (
    pd.DataFrame(combined_df.target.str.split("|").tolist(), index=combined_df.index)
    .stack()
    .reset_index()
)

target_split_df.columns = [split_col_index, "_", "target_unique"]

print(target_split_df.shape)
target_split_df.head()


# In[18]:


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
combined_df.to_csv(output_file, sep='\t', index=False)

print(long_combined_df.shape)
long_combined_df.head()

