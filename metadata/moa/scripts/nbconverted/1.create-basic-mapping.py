#!/usr/bin/env python
# coding: utf-8

# # Create a Basic Mapping of Broad IDs to MOAs and Targets
# 
# Here, I add the `pert_id` column to an additional file.
# This file contains only the essential mapping columns between cell painting data and drug annotations.
# 
# The primary addition is truncating the `broad_id` column into the `pert_id` column (e.g. `BRD-K89787693-001-01-1` becomes `BRD-K89787693`) and describing the impact.
# 
# I create a new file (called `repurposing_info_basic.tsv`) that only contains the unique columns `pert_id`, `pert_iname`, `moa`, and `target`.
# 
# The `broad_id` column contains additional supplier, batch and aliquot information.
# The `pert_id` column is the essential information that is used to directly map compounds to profiles.

# In[1]:


import os
import pandas as pd


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


use_cols = ["pert_id", "pert_iname", "moa", "target"]


# ## Load Data and Create the `pert_id` column

# In[4]:


annotation_df = pd.read_csv("repurposing_info.tsv", sep='\t')

annotation_df = annotation_df.assign(
    pert_id=annotation_df.broad_id.str.slice(0, 13)
)

print(annotation_df.shape)
annotation_df.head()


# ## Describe the Effect of adding `pert_id`
# 
# How many columns are replicated?

# In[5]:


# There are no duplicate `broad_ids`
annotation_df.broad_id.duplicated().sum()


# In[6]:


pert_counts = annotation_df.pert_id.value_counts()
pert_counts.head(10)


# In[7]:


pert_counts.hist();


# In[8]:


# How many pert ids are duplicated?
num_duplicate = (pert_counts > 1).sum()
percent_duplicate = (num_duplicate / pert_counts.shape[0]) * 100

print("There are {} ({}%) duplicated `pert_ids`".format(num_duplicate, percent_duplicate))


# In[9]:


# Examples of duplicate pert_id columns
top_duplicated_pert_id = pert_counts.head(1).index.values[0]
annotation_df.query("pert_id == @top_duplicated_pert_id")


# ## Select only essential columns and determine/reconcile any discrepancies

# In[10]:


basic_df = annotation_df.loc[:, use_cols].drop_duplicates()

print(basic_df.shape)
basic_df.head()


# In[11]:


duplicated_pert_ids = basic_df.pert_id.loc[basic_df.pert_id.duplicated()]
assert len(duplicated_pert_ids) == 0, "Warning! There are duplicated pert_ids"

