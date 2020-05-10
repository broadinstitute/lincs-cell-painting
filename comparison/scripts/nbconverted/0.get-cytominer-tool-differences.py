#!/usr/bin/env python
# coding: utf-8

# # Comparing Pycytominer and Cytominer Processing
# 
# We have previously processed all of the Drug Repurposing Hub Cell Painting Data using [cytominer](https://github.com/cytomining/cytominer).
# Cytominer is an R based image-based profiling tool.
# In this repo, we reprocess the data with [pycytominer](https://github.com/cytomining/pycytominer).
# As the name connotes, pycytominer is a python based image-based profiling tool.
# 
# We include all processing scripts and present the pycytominer profiles in this open source repository.
# The repository represents a unified bioinformatics pipeline applied to all Cell Painting Drug Repurposing Profiles. In this notebook, we compare the resulting output data between the processing pipelines for the two tools: Cytominer and pycytominer.
# We output several metrics comparing the two approaches
# 
# ## Metrics
# 
# In all cases, we calculate the element-wise absolute value difference between pycytominer and cytominer profiles.
# 
# 1. Mean, median, and sum of element-wise differencs
# 2. Per feature mean, median, and sum of element-wise differences
# 3. Feature selection procedure differences per feature (level 4b only)
# 
# In addition, we confirm alignment of the following metadata columns:
# 
# * Well
# * Broad Sample Name
# * Plate
# 
# Other metadata columns are not expected to be aligned.
# For example, we have [updated MOA and Target information](https://github.com/broadinstitute/lincs-cell-painting/issues/11) in the pycytominer version.
# 
# ## Data Levels
# 
# Image-based profiling results in the following output data levels.
# We do not compare all data levels in this notebook.
# 
# | Data | Level | Comparison |
# | :---- | :---- | :-------- |
# | Images | Level 1 | NA |
# | SQLite File (single cell profiles ) | Level 2 | NA |
# | Aggregated Profiles with Well Information (metadata) | Level 3 | Yes  |
# | Normalized Aggregated Profiles with Metadata | Level 4a | Yes | 
# | Normalized and Feature Selected Aggregated Profiles with Metadata | Level 4b | Yes |
# | Perturbation Profiles created Summarizing Replicates  | Level 5 | No |

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd

from util import build_file_dictionary, load_data, generate_output_filenames


# In[3]:


# Pycytominer plates are saved with 5 floating point decimals
round_decimals = 5

# Create the output directory
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)


# In[4]:


def get_metrics(pycyto_df, cyto_df, features):
    # Align features
    pycyto_df = pycyto_df.reindex(features, axis="columns")
    cyto_df = cyto_df.reindex(features, axis="columns")

    # Assess difference
    abs_diff = pycyto_df.subtract(cyto_df).abs()
    mean_diff = abs_diff.mean()
    median_diff = abs_diff.median()
    sum_diff = abs_diff.sum()

    complete_mean_diff = mean_diff.replace([np.inf, -np.inf], np.nan).dropna().mean()
    complete_median_diff = (
        median_diff.replace([np.inf, -np.inf], np.nan).dropna().mean()
    )
    complete_sum_diff = sum_diff.replace([np.inf, -np.inf], np.nan).dropna().sum()

    return (
        mean_diff,
        complete_mean_diff,
        median_diff,
        complete_median_diff,
        sum_diff,
        complete_sum_diff,
    )


def find_feature_diff(pycyto_df, cyto_df, plate, all_features):
    all_features_df = pd.DataFrame(
        ["missing"] * len(all_features), index=all_features, columns=[plate]
    )
    pycyto_features = set(pycyto_df.columns.tolist())
    cyto_features = set(cyto_df.columns.tolist())
    present_both = pycyto_features.intersection(cyto_features)

    all_features_df.loc[
        all_features_df.index.isin(pycyto_features), plate
    ] = "only_pycytominer"
    all_features_df.loc[
        all_features_df.index.isin(cyto_features), plate
    ] = "only_cytominer"
    all_features_df.loc[
        all_features_df.index.isin(present_both), plate
    ] = "present_both"

    return all_features_df


# In[5]:


pycytominer_dir = pathlib.Path("../profiles/backend/")
cytominer_dir = pathlib.Path(
    "/Users/gway/work/projects/"
    + "2015_10_05_DrugRepurposing_AravindSubramanian_GolubLab_Broad"
    + "/workspace/backend/2016_04_01_a549_48hr_batch1/"
)


# In[6]:


pycytominer_plate_files = build_file_dictionary(pycytominer_dir, tool="pycytominer")
cytominer_plate_files = build_file_dictionary(cytominer_dir, tool="cytominer")


# In[7]:


pycytominer_plates = set(sorted(pycytominer_plate_files.keys()))
cytominer_plates = set(sorted(cytominer_plate_files.keys()))

assert (
    cytominer_plates == pycytominer_plates
), "Stop, not every plate is measured using both tools"

print(len(pycytominer_plates))


# In[8]:


level_3_mean_diff = []
level_3_completemean_diff = {}
level_3_median_diff = []
level_3_completemedian_diff = {}
level_3_sum_diff = []
level_3_completesum_diff = {}

level_4a_mean_diff = []
level_4a_completemean_diff = {}
level_4a_median_diff = []
level_4a_completemedian_diff = {}
level_4a_sum_diff = []
level_4a_completesum_diff = {}

level_4b_mean_diff = []
level_4b_completemean_diff = {}
level_4b_sum_diff = []
level_4b_median_diff = []
level_4b_completemedian_diff = {}
level_4b_completesum_diff = {}
level_4b_feature_select = []

# Calculate metrics per plate
for plate in list(cytominer_plates):
    # Calculate level 3 metrics
    pycyto_df, cyto_df = load_data(
        plate,
        pycytominer_plate_files,
        cytominer_plate_files,
        level="level_3",
        round_decimals=round_decimals,
    )
    # Define features (note that the features were checked and aligned in load_data)
    features = pycyto_df.columns.tolist()
    # Get differences
    (
        mean_diff,
        complete_mean_diff,
        median_diff,
        complete_median_diff,
        sum_diff,
        complete_sum_diff,
    ) = get_metrics(pycyto_df, cyto_df, features)
    # Store results
    level_3_mean_diff.append(mean_diff)
    level_3_completemean_diff[plate] = complete_mean_diff
    level_3_median_diff.append(median_diff)
    level_3_completemedian_diff[plate] = complete_median_diff
    level_3_sum_diff.append(sum_diff)
    level_3_completesum_diff[plate] = complete_sum_diff

    # Calculate level 4a metrics
    pycyto_df, cyto_df = load_data(
        plate,
        pycytominer_plate_files,
        cytominer_plate_files,
        level="level_4a",
        round_decimals=round_decimals,
    )
    # Get differences
    (
        mean_diff,
        complete_mean_diff,
        median_diff,
        complete_median_diff,
        sum_diff,
        complete_sum_diff,
    ) = get_metrics(pycyto_df, cyto_df, features)
    # Store results
    level_4a_mean_diff.append(mean_diff)
    level_4a_completemean_diff[plate] = complete_mean_diff
    level_4a_median_diff.append(median_diff)
    level_4a_completemedian_diff[plate] = complete_median_diff
    level_4a_sum_diff.append(sum_diff)
    level_4a_completesum_diff[plate] = complete_sum_diff

    # Calculate level 4b metrics
    pycyto_df, cyto_df = load_data(
        plate,
        pycytominer_plate_files,
        cytominer_plate_files,
        level="level_4b",
        round_decimals=round_decimals,
    )
    # Determine feature selection differences
    feature_select_df = find_feature_diff(pycyto_df, cyto_df, plate, features)
    features_present_in_both = feature_select_df.loc[
        feature_select_df.loc[:, plate] == "present_both", plate
    ].index.tolist()
    # Get differences
    (
        mean_diff,
        complete_mean_diff,
        median_diff,
        complete_median_diff,
        sum_diff,
        complete_sum_diff,
    ) = get_metrics(pycyto_df, cyto_df, features_present_in_both)
    # Store results
    level_4b_mean_diff.append(mean_diff)
    level_4b_completemean_diff[plate] = complete_mean_diff
    level_4b_median_diff.append(median_diff)
    level_4b_completemedian_diff[plate] = complete_median_diff
    level_4b_sum_diff.append(sum_diff)
    level_4b_completesum_diff[plate] = complete_sum_diff
    level_4b_feature_select.append(feature_select_df)


# ## Compile Results

# In[9]:


level_3_mean_diff_df = pd.concat(level_3_mean_diff, axis="columns", sort=True)
level_3_mean_diff_df.columns = list(cytominer_plates)
level_3_completemean_diff_df = pd.DataFrame(
    level_3_completemean_diff, index=["complete_mean_diff"]
).transpose()

level_3_median_diff_df = pd.concat(level_3_median_diff, axis="columns", sort=True)
level_3_median_diff_df.columns = list(cytominer_plates)
level_3_completemedian_diff_df = pd.DataFrame(
    level_3_completemedian_diff, index=["complete_median_diff"]
).transpose()

level_3_sum_diff_df = pd.concat(level_3_sum_diff, axis="columns", sort=True)
level_3_sum_diff_df.columns = list(cytominer_plates)
level_3_completesum_diff_df = pd.DataFrame(
    level_3_completesum_diff, index=["complete_sum_diff"]
).transpose()


# In[10]:


level_4a_mean_diff_df = pd.concat(level_4a_mean_diff, axis="columns")
level_4a_mean_diff_df.columns = list(cytominer_plates)
level_4a_completemean_diff_df = pd.DataFrame(
    level_4a_completemean_diff, index=["complete_mean_diff"]
).transpose()

level_4a_median_diff_df = pd.concat(level_4a_median_diff, axis="columns", sort=True)
level_4a_median_diff_df.columns = list(cytominer_plates)
level_4a_completemedian_diff_df = pd.DataFrame(
    level_4a_completemedian_diff, index=["complete_median_diff"]
).transpose()

level_4a_sum_diff_df = pd.concat(level_4a_sum_diff, axis="columns", sort=True)
level_4a_sum_diff_df.columns = list(cytominer_plates)
level_4a_completesum_diff_df = pd.DataFrame(
    level_4a_completesum_diff, index=["complete_sum_diff"]
).transpose()


# In[11]:


level_4b_mean_diff_df = pd.concat(level_4b_mean_diff, axis="columns")
level_4b_mean_diff_df.columns = list(cytominer_plates)
level_4b_completemean_diff_df = pd.DataFrame(
    level_4b_completemean_diff, index=["complete_mean_diff"]
).transpose()

level_4b_median_diff_df = pd.concat(level_4b_median_diff, axis="columns", sort=True)
level_4b_median_diff_df.columns = list(cytominer_plates)
level_4b_completemedian_diff_df = pd.DataFrame(
    level_4b_completemedian_diff, index=["complete_median_diff"]
).transpose()

level_4b_sum_diff_df = pd.concat(level_4b_sum_diff, axis="columns", sort=True)
level_4b_sum_diff_df.columns = list(cytominer_plates)
level_4b_completesum_diff_df = pd.DataFrame(
    level_4b_completesum_diff, index=["complete_sum_diff"]
).transpose()

level_4b_feature_select_df = pd.concat(level_4b_feature_select, axis="columns")


# ## Output Results

# In[12]:


level = "level_3"
level_3_files = generate_output_filenames(output_dir, level)

# Output mean
level_3_mean_diff_df.to_csv(
    level_3_files["mean"], sep="\t", index=True, compression="gzip"
)

# Output median
level_3_median_diff_df.to_csv(
    level_3_files["median"], sep="\t", index=True, compression="gzip"
)

# Output sum
level_3_sum_diff_df.to_csv(
    level_3_files["sum"], sep="\t", index=True, compression="gzip"
)


# In[13]:


level = "level_4a"
level_4a_files = generate_output_filenames(output_dir, level)

# Output mean
level_4a_mean_diff_df.to_csv(
    level_4a_files["mean"], sep="\t", index=True, compression="gzip"
)

# Output median
level_4a_median_diff_df.to_csv(
    level_4a_files["median"], sep="\t", index=True, compression="gzip"
)

# Output sum
level_4a_sum_diff_df.to_csv(
    level_4a_files["sum"], sep="\t", index=True, compression="gzip"
)


# In[14]:


level = "level_4b"
level_4b_files = generate_output_filenames(output_dir, level)

# Output mean
level_4b_mean_diff_df.to_csv(
    level_4b_files["mean"], sep="\t", index=True, compression="gzip"
)

# Output median
level_4b_median_diff_df.to_csv(
    level_4b_files["median"], sep="\t", index=True, compression="gzip"
)

# Output sum
level_4b_sum_diff_df.to_csv(
    level_4b_files["sum"], sep="\t", index=True, compression="gzip"
)

# Output feature select summary file
output_file = f"{output_dir}/comparison_result_4b_feature_select.tsv.gz"
level_4b_feature_select_df.to_csv(output_file, sep="\t", index=True, compression="gzip")


# In[15]:


# Concatenate level 3 results
level_3_complete_df = pd.concat(
    [
        level_3_completemean_diff_df,
        level_3_completemedian_diff_df,
        level_3_completesum_diff_df,
    ],
    axis="columns",
)

level_3_complete_df.columns = [f"level_3_{x}" for x in level_3_complete_df.columns]

# Concatenate level 4a results
level_4a_complete_df = pd.concat(
    [
        level_4a_completemean_diff_df,
        level_4a_completemedian_diff_df,
        level_4a_completesum_diff_df,
    ],
    axis="columns",
)

level_4a_complete_df.columns = [f"level_4a_{x}" for x in level_4a_complete_df.columns]

# Concatenate level 4b results
level_4b_complete_df = pd.concat(
    [
        level_4b_completemean_diff_df,
        level_4b_completemedian_diff_df,
        level_4b_completesum_diff_df,
    ],
    axis="columns",
)

level_4b_complete_df.columns = [f"level_4b_{x}" for x in level_4b_complete_df.columns]

# Combine all results
complete_df = pd.concat(
    [level_3_complete_df, level_4a_complete_df, level_4b_complete_df,], axis="columns",
)

# Output file
output_file = f"{output_dir}/comparison_result_metric_summary.tsv"
complete_df.to_csv(output_file, sep="\t", index=True)

print(complete_df.shape)
complete_df.head()

