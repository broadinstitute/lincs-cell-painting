#!/usr/bin/env python
# coding: utf-8

# # Align Platemap and MOA Data
# 
# Here, I align Broad IDs in the Platemap Data to the Most Recent MOA/Target Info

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[3]:


def split_pipe_broad_id(map_subset_df, broad_id_col):
    # Split a dataframe with broad_id_cols that have pipe delimiters
    broad_split_df = (
        pd.DataFrame(
            map_subset_df.loc[:, broad_id_col]
            .squeeze()
            .dropna()
            .str.split("|")
            .tolist(),
            index=map_df.dropna(subset=[broad_id_col]).index,
        )
        .stack()
        .reset_index()
    )
    broad_split_df.columns = ["map_index", "_", "broad_id_unique"]

    map_long_df = (
        map_subset_df.merge(
            broad_split_df.loc[:, ["map_index", "broad_id_unique"]],
            left_index=True,
            right_on="map_index",
            how="left",
        )
        .reset_index(drop=True)
        .loc[:, ["InChIKey14", "broad_id_unique"]]
        .drop_duplicates()
        .assign(broad_date=broad_id_col)
        .dropna()
        .reset_index(drop=True)
    )

    return map_long_df


# ## Step 1: Define InChIKey14-to-Broad ID Map

# In[4]:


mapping_col = "InChIKey14"

repurposing_col = "broad_id_20170327"
broad_id_other_dates = ["20170327", "20180516", "20180907", "20200324"]


# In[5]:


# Load platemap data
file = pathlib.PurePath("../platemaps/broad_sample_info.tsv")
all_platemap_df = pd.read_csv(file, sep="\t").dropna(
    subset=["broad_sample", "broad_id"]
)

print(all_platemap_df.shape)
all_platemap_df.head()


# In[6]:


# How many unique broad IDs?
len(all_platemap_df.broad_id.unique())


# In[7]:


# Load MOA map
map_file = pathlib.PurePath("clue/broad_id_map.csv")
map_df = pd.read_csv(map_file)

print(map_df.shape)
map_df.head()


# In[8]:


# Create MOA map long mapping every Broad ID to InChIKey14
all_broad_map = []
for date in broad_id_other_dates:
    broad_id_col = f"broad_id_{date}"
    deprecated_broad_id_col = f"deprecated_broad_id_{date}"

    select_cols = [mapping_col, broad_id_col, deprecated_broad_id_col]
    map_subset_df = map_df.loc[:, select_cols]

    broad_long_df = split_pipe_broad_id(map_subset_df, broad_id_col)
    all_broad_map.append(broad_long_df)

    if date != "20170327":
        deprecated_broad_long_df = split_pipe_broad_id(
            map_subset_df, deprecated_broad_id_col
        )
        all_broad_map.append(deprecated_broad_long_df)

# Combine full broad_id map
all_broad_map = (
    pd.concat(all_broad_map)
    # This step drops InChIKeys and Broad_ID combos that exist in multiple dates
    .drop_duplicates(subset=[mapping_col, "broad_id_unique"]).reset_index(drop=True)
)

print(all_broad_map.shape)
all_broad_map.head()


# In[9]:


# Are there Duplicate InChIKey14s?
inchi_counts = all_broad_map.loc[:, mapping_col].squeeze().value_counts()
inchi_counts.hist()


# In[10]:


# Show an example of the most egregious multi-mapping InChiKey14
all_broad_map.loc[all_broad_map.loc[:, mapping_col] == inchi_counts.index[0], :]


# ## Step 2: Perform the Merge with the Platemap Data
# 
# This step identifies the InChIKey14 map for the Cell Painting Repurposing Hub Data Subset.

# In[11]:


merged_platemap_df = (
    all_platemap_df.merge(
        all_broad_map, left_on="broad_id", right_on="broad_id_unique", how="left"
    )
    # We don't care about plate info here
    .drop_duplicates(subset=["broad_id", mapping_col]).loc[
        :, ["broad_sample", "broad_id", mapping_col, "broad_date"]
    ]
    # Drop DMSO
    .dropna(subset=["broad_sample"])
)

print(merged_platemap_df.shape)
merged_platemap_df.head()


# In[12]:


# Are there Duplicate Broad IDs
broad_id_counts = merged_platemap_df.broad_id.value_counts()
broad_id_counts.hist()


# In[13]:


# There are indeed duplicated broad IDs!
# This means that some of these broad IDs will recieve alternative MOA/Targets
duplicate_broad_ids = broad_id_counts[broad_id_counts > 1].index.tolist()
merged_platemap_df.query("broad_id in @duplicate_broad_ids")


# ## Step 3: Align Repurposing Platemap with Most Recent MOA/Target Info

# In[14]:


# Load Data
moa_file = "repurposing_info.tsv"
moa_df = pd.read_csv(moa_file, sep="\t")

print(moa_df.shape)
moa_df.head()


# In[15]:


# Perform the alignment
platemap_moa_df = merged_platemap_df.merge(
    moa_df.loc[:, [mapping_col, "pert_iname", "moa", "target"]].drop_duplicates(),
    on=mapping_col,
    how="left",
).loc[:, ["broad_sample", "broad_id", mapping_col, "moa", "target", "broad_date"]]

print(platemap_moa_df.shape)
platemap_moa_df.head()


# In[16]:


# Show an example of a broad_id that must be resolved
platemap_moa_df.query("broad_id == 'BRD-K56844688'")


# ## Step 5 - Quantify what's missing

# In[17]:


# Some samples are missing an InChIKey14
platemap_moa_with_inchi_df = platemap_moa_df.dropna(subset=[mapping_col])

print(platemap_moa_with_inchi_df.shape)
platemap_moa_with_inchi_df.head()


# In[18]:


missing_inchi_key = set(platemap_moa_df.broad_id.unique()).difference(
    set(platemap_moa_with_inchi_df.broad_id)
)
len(missing_inchi_key)


# In[19]:


all_unique_broad_platemap_ids = set(platemap_moa_df.broad_id.unique())
len(all_unique_broad_platemap_ids)


# In[20]:


complete_info = set(
    platemap_moa_with_inchi_df.loc[
        ~(
            platemap_moa_with_inchi_df.moa.isna()
            | platemap_moa_with_inchi_df.target.isna()
        ),
        "broad_id",
    ].values
)
len(complete_info)


# In[21]:


moa_missing = set(
    platemap_moa_with_inchi_df.sort_values(by="moa", na_position="last")
    .drop_duplicates(subset="broad_id", keep="first")
    .loc[platemap_moa_with_inchi_df.moa.isna(), "broad_id"]
    .values
).difference(complete_info)

moa_present = all_unique_broad_platemap_ids.difference(moa_missing)
len(moa_missing)


# In[22]:


target_missing = set(
    platemap_moa_with_inchi_df.sort_values(by="target", na_position="last")
    .drop_duplicates(subset="broad_id", keep="first")
    .loc[platemap_moa_with_inchi_df.target.isna(), "broad_id"]
    .values
).difference(complete_info)

target_present = all_unique_broad_platemap_ids.difference(target_missing)
len(target_missing)


# In[23]:


missing_moa_info = moa_missing.difference(target_missing)
len(missing_moa_info)


# In[24]:


missing_target_info = target_missing.difference(moa_missing)
len(missing_target_info)


# In[25]:


missing_both_info = moa_missing.intersection(target_missing)

missing_both_info = (
    missing_both_info.difference(complete_info)
    .difference(target_present)
    .difference(moa_present)
)
len(missing_both_info)


# In[26]:


pd.DataFrame(
    {
        "Unique Broad IDs": len(all_unique_broad_platemap_ids),
        "Broad IDs with Complete MOA/Target": len(complete_info),
        "Missing InChIKey14": len(missing_inchi_key),
        "Broad IDs with Only Missing MOA Info": len(missing_moa_info),
        "Broad IDs with Only Missing Target Info": len(missing_target_info),
        "Broad IDs Missing Both MOA and Target": len(missing_both_info),
    },
    index=["count"],
).transpose()


# ## Step 6 - Resolve Stereochemistry Differences

# In[27]:


alternative_delim = ";"


# In[28]:


# Drop rows in the following scenario:
# 1. The broad_id is assigned with two or more InChIKey14
# 2. Drop InChIKey14s in these rows that have missing values in both moa and target

broad_id_with_multiple_inchi = (
    platemap_moa_with_inchi_df.groupby("broad_id")[mapping_col]
    .count()
    .sort_values(ascending=False)
)

broad_id_with_singular_inchi = set(
    broad_id_with_multiple_inchi[broad_id_with_multiple_inchi == 1].index
)
print(f"Only one InChIKey14: {len(broad_id_with_singular_inchi)}")

broad_id_with_multiple_inchi = set(
    broad_id_with_multiple_inchi[broad_id_with_multiple_inchi > 1].index
)
print(f"More than one InChIKey14: {len(broad_id_with_multiple_inchi)}")


# In[29]:


# First, drop entries that have missing values in BOTH moa and target
# Also drop duplicated columns
platemap_multiple_inchi_select_df = platemap_moa_with_inchi_df.query(
    "broad_id in @broad_id_with_multiple_inchi"
)

platemap_multiple_inchi_select_df = platemap_multiple_inchi_select_df.loc[
    ~(
        platemap_multiple_inchi_select_df.moa.isna()
        & platemap_multiple_inchi_select_df.target.isna()
    )
].drop_duplicates()

print(platemap_multiple_inchi_select_df.shape)
platemap_multiple_inchi_select_df.head(20)

assert len(broad_id_with_multiple_inchi) == len(
    platemap_multiple_inchi_select_df.broad_id.unique()
), "Stop, a broad ID was dropped somewhere that shouldn't have"


# In[30]:


# After this step, are there any resolved broad_ids?
first_step_resolve = platemap_multiple_inchi_select_df.broad_id.value_counts()
first_step_resolve = list(first_step_resolve[first_step_resolve == 1].index)

first_step_resolve_df = platemap_multiple_inchi_select_df.query(
    "broad_id in @first_step_resolve"
).reset_index(drop=True)

print(first_step_resolve_df.shape)
first_step_resolve_df.head()


# In[31]:


# Remove these from other broad IDs that need to be resolved
platemap_alternative_inchi_df = platemap_multiple_inchi_select_df.query(
    "broad_id not in @first_step_resolve"
).drop_duplicates(["broad_id", "moa", "target"])

platemap_alternative_inchi_df = platemap_alternative_inchi_df.assign(
    deprecated=[x[0] for x in platemap_alternative_inchi_df.broad_date.str.split("_")]
)

assert (
    len(platemap_alternative_inchi_df.deprecated.unique()) == 1
), "Warning! Deprecated info should not be resolved"

print(platemap_alternative_inchi_df.shape)
platemap_alternative_inchi_df.head()


# In[32]:


# Now resolve the remaining broad_ids
alternative_resolved_df = []
for broad_id in platemap_alternative_inchi_df.broad_id.unique():
    subset_platemap_alternative_inchi_df = platemap_alternative_inchi_df.query(
        "broad_id == @broad_id"
    ).reset_index(drop=True)

    primary_moa_df = subset_platemap_alternative_inchi_df.iloc[0, :]
    alternative_moa_df = subset_platemap_alternative_inchi_df.drop(
        primary_moa_df.name, axis="rows"
    )

    moa_join = alternative_delim.join([str(x) for x in alternative_moa_df.moa.tolist()])
    target_join = alternative_delim.join(
        [str(x) for x in alternative_moa_df.target.tolist()]
    )
    alternative_resolved_df.append(
        pd.DataFrame(primary_moa_df)
        .transpose()
        .assign(alternative_moa=moa_join, alternative_target=target_join)
    )

alternative_resolved_df = pd.concat(alternative_resolved_df)

assert sorted(alternative_resolved_df.broad_id) == sorted(
    platemap_alternative_inchi_df.broad_id.unique()
), "Warning, something went missing along the way"

assert (
    alternative_resolved_df.broad_id.value_counts().sum()
    == alternative_resolved_df.shape[0]
), "Warning, there are still duplicate broad_ids"

print(alternative_resolved_df.shape)
alternative_resolved_df.head()


# ## Step 7 - Combine Annotated Platemap with Resolution

# In[33]:


# First, remove the broad_ids that needed to be resolved
resolved_broad_ids = (
    alternative_resolved_df.broad_id.tolist() + first_step_resolve_df.broad_id.tolist()
)

platemap_moa_easy_df = platemap_moa_df.query("broad_id not in @resolved_broad_ids")

assert (
    platemap_moa_easy_df.broad_id.value_counts().sum() == platemap_moa_easy_df.shape[0]
), "Error, there are still broad_ids that require resolution"


# In[34]:


platemap_moa_easy_df = platemap_moa_easy_df.assign(
    alternative_moa=np.nan, alternative_target=np.nan
)

first_step_resolve_df = first_step_resolve_df.assign(
    alternative_moa=np.nan, alternative_target=np.nan
)


# In[35]:


moa_map_df = (
    pd.concat([platemap_moa_easy_df, first_step_resolve_df, alternative_resolved_df])
    .sort_values(by="broad_id")
    .reset_index(drop=True)
)

assert (
    moa_map_df.broad_id.value_counts().sum() == moa_map_df.shape[0]
), "Error, broad_ids have been duplicated somewhere along the way"

print(moa_map_df.shape)
moa_map_df.head()


# In[36]:


# Interpretation: This is the CLUE version the 2020 CLUE version (most recent) used to map
moa_map_df.broad_date.value_counts()


# In[37]:


map_output_file = "repurposing_info_external_moa_map_resolved.tsv"
moa_map_df.to_csv(map_output_file, sep="\t", index=False)


# In[38]:


pd.DataFrame(
    {
        "All Unique Broad IDs": len(all_unique_broad_platemap_ids),
        "One-to-One Mapping": platemap_moa_easy_df.query(
            "broad_id not in list(@missing_inchi_key)"
        ).shape[0],
        "Missing InChiKey14": len(missing_inchi_key),
        "Duplicate InChIKey14 with both MOA/Target Missing": first_step_resolve_df.shape[
            0
        ],
        "Duplicate InChIKey14 requires alternative MOA/Target annotations": alternative_resolved_df.shape[
            0
        ],
    },
    index=["count"],
).transpose()

