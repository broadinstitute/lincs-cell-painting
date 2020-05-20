#!/usr/bin/env python
# coding: utf-8

# # Summarize Tool Differences
# 
# In the following notebook, we summarize the differences we calculated in `0.get-cytominer-tool-differences`.
# We summarize the results in a series of visualizations and descriptive statistics.

# In[1]:


get_ipython().run_line_magic('load_ext', 'nb_black')


# In[2]:


import os
import pathlib
import numpy as np
import pandas as pd
import plotnine as gg

from util import build_filenames


# In[3]:


# Set constants
batch = "2016_04_01_a549_48hr_batch1"
input_dir = pathlib.Path("results", batch)
levels = ["level_3", "level_4a", "level_4b", "pycytominer_select"]
metrics = ["mean", "median", "sum"]

# Set output directory
output_fig_dir = pathlib.Path("figures", batch)
output_fig_dir.mkdir(parents=True, exist_ok=True)

# Set plotting defaults
dpi = 500
height = 3.5
width = 6

# Set common plotnine theme
theme_summary = gg.theme_bw() + gg.theme(
    axis_text_x=gg.element_blank(),
    axis_text_y=gg.element_text(size=6),
    axis_title=gg.element_text(size=8),
    strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
)


# In[4]:


# Load Data
results_files = {}
for level in levels:
    file_names = build_filenames(input_dir, level=level)
    metric_df = {}
    for metric in file_names:
        df = pd.read_csv(file_names[metric], sep="\t", index_col=0)
        metric_df[metric] = df

    results_files[level] = metric_df

summary_file = pathlib.Path(f"{input_dir}/comparison_result_metric_summary.tsv")
summary_df = (
    pd.read_csv(summary_file, sep="\t", index_col=0)
    .sort_index()
    .reset_index()
    .rename({"index": "plate"}, axis="columns")
)

print(summary_df.shape)
summary_df.head()


# In[5]:


# Isolate the outlier plates that were processed by Cytominer in a different way
# See https://github.com/broadinstitute/lincs-cell-painting/issues/3#issuecomment-591994451
nonuniform_plates = summary_df.query("level_3_complete_median_diff > 1").plate.tolist()
nonuniform_plates


# In[6]:


# Since we know that these were processed differently, remove from the comparison
summary_df = summary_df.query("plate not in @nonuniform_plates").reset_index(drop=True)


# In[7]:


# Ensure the plates are in order in each plot
plate_order = summary_df.plate.tolist()
summary_df.plate = pd.Categorical(
    summary_df.plate, categories=plate_order, ordered=True
)


# In[8]:


# Process summary dataframe
summary_melted_df = []
for level in levels:
    # Summary for pycytominer select not shown to reduce comparison size
    if level == "pycytominer_select":
        continue
    for metric in metrics:
        col_name = f"{level}_complete_{metric}_diff"
        subset_df = (
            summary_df.loc[:, ["plate", col_name]]
            .assign(metric=metric, level=level)
            .rename({col_name: "metric_value"}, axis="columns")
        )
        summary_melted_df.append(subset_df)

summary_melted_df = pd.concat(summary_melted_df).reset_index(drop=True)
summary_melted_df = summary_melted_df.assign(uniform=True)
summary_melted_df.loc[
    summary_melted_df.plate.isin(nonuniform_plates), "uniform"
] = False

print(summary_melted_df.shape)
summary_melted_df.head()


# In[9]:


summary_metric_full_gg = (
    gg.ggplot(summary_melted_df, gg.aes(x="plate", y="metric_value"))
    + gg.geom_point(size=0.5)
    + gg.facet_grid("metric~level", scales="free")
    + gg.xlab("Plate")
    + gg.ylab("Absolute Difference\nBetween Tools")
    + gg.ggtitle("Per Plate Summary (Full)")
    + theme_summary
)

output_file = pathlib.Path(f"{output_fig_dir}/summary_metrics_full.png")
summary_metric_full_gg.save(output_file, dpi=500, height=4, width=6)

print(summary_metric_full_gg)

summary_metric_zoom_gg = (
    gg.ggplot(
        summary_melted_df.query("metric_value < 10"),
        gg.aes(x="plate", y="metric_value"),
    )
    + gg.geom_point(size=0.5)
    + gg.facet_grid("metric~level", scales="free")
    + gg.xlab("Plate")
    + gg.ylab("Absolute Difference\nBetween Tools")
    + gg.ggtitle("Per Plate Summary (Zoom)")
    + theme_summary
)

output_file = pathlib.Path(f"{output_fig_dir}/summary_metrics_zoom.png")
summary_metric_zoom_gg.save(output_file, dpi=500, height=3.5, width=6)

print(summary_metric_zoom_gg)


# In[10]:


# Wrangle output metric data to be plot ready
all_feature_results_df = []
for level in levels:
    for metric in metrics:
        plot_ready_df = (
            results_files[level][metric]
            .reset_index()
            .rename({"index": "feature"}, axis="columns")
            .melt(id_vars="feature", var_name="plate", value_name="metric_value")
            .assign(metric=metric, level=level)
        )
        all_feature_results_df.append(plot_ready_df)

all_feature_results_df = pd.concat(all_feature_results_df).reset_index(drop=True)

# Predetermine feature order
feature_order = sorted(list(set(all_feature_results_df.feature)))

all_feature_results_df = all_feature_results_df.assign(uniform=True)
all_feature_results_df.loc[
    all_feature_results_df.plate.isin(nonuniform_plates), "uniform"
] = False

all_feature_results_df.plate = pd.Categorical(
    all_feature_results_df.plate, categories=plate_order, ordered=True
)
all_feature_results_df.feature = pd.Categorical(
    all_feature_results_df.feature, categories=feature_order, ordered=True
)

# Select only the uniform plates
all_feature_results_df = all_feature_results_df.query("uniform")

print(all_feature_results_df.shape)
all_feature_results_df.head()


# ## All Feature and Plate Summary

# In[11]:


all_feature_results_df.groupby(["metric", "level"])["metric_value"].describe()


# ## Generate Two Figures Per Data Level
# 
# Split into different cells to prevent kernel death.

# In[12]:


for level in levels:
    all_feature_results_subset_df = all_feature_results_df.query(
        "level == @level"
    ).reset_index(drop=True)

    output_dir = pathlib.Path(f"{output_fig_dir}/{level}")
    output_dir.mkdir(exist_ok=True)

    # Figure 1 - Per plate feature differences
    per_plate_feature_gg = (
        gg.ggplot(all_feature_results_subset_df, gg.aes(x="plate", y="metric_value"))
        + gg.geom_point(size=0.1, alpha=0.5)
        + gg.facet_wrap("~metric", scales="free", nrow=len(metrics))
        + gg.xlab("Plate")
        + gg.ylab("Feature Difference\nBetween Tools")
        + gg.ggtitle(f"Plate Summary\n{level}")
        + theme_summary
    )

    output_file = pathlib.Path(f"{output_dir}/{level}_metrics_per_plate.png")
    per_plate_feature_gg.save(output_file, dpi=dpi, height=height, width=width)

    print(per_plate_feature_gg)
    del per_plate_feature_gg


# In[13]:


for level in levels:
    all_feature_results_subset_df = all_feature_results_df.query(
        "level == @level"
    ).reset_index(drop=True)

    output_dir = pathlib.Path(f"{output_fig_dir}/{level}")
    output_dir.mkdir(exist_ok=True)

    # Figure 2 - Per feature plate differences
    per_feature_gg = (
        gg.ggplot(all_feature_results_subset_df, gg.aes(x="feature", y="metric_value"))
        + gg.geom_point(size=0.1, alpha=0.5)
        + gg.facet_wrap("~metric", scales="free", nrow=len(metrics))
        + gg.xlab("Feature")
        + gg.ylab("Feature Difference\nBetween Tools")
        + gg.ggtitle(f"Feature Summary Across Plates\n{level}")
        + theme_summary
    )

    output_file = pathlib.Path(f"{output_dir}/{level}_metrics_per_feature.png")
    per_feature_gg.save(output_file, dpi=dpi, height=height, width=width)

    print(per_feature_gg)
    del per_feature_gg


# ## Feature Selection Summary

# In[14]:


# Load data
select_file = pathlib.Path(f"{input_dir}/comparison_result_4b_feature_select.tsv.gz")
select_df = (
    pd.read_csv(select_file, sep="\t", index_col=0)
    .reset_index()
    .rename({"index": "feature"}, axis="columns")
    .melt(id_vars="feature", var_name="plate", value_name="status")
    .query("plate not in @nonuniform_plates")
)

# Reorder data
select_df.plate = pd.Categorical(select_df.plate, categories=plate_order, ordered=True)
select_df.feature = pd.Categorical(
    select_df.feature, categories=feature_order, ordered=True
)

print(select_df.shape)
select_df.head()


# In[15]:


feature_select_gg = (
    gg.ggplot(select_df, gg.aes(x="feature", y="plate", fill="status"))
    + gg.geom_tile(size=0.5)
    + gg.ggtitle("Feature Select Summary")
    + theme_summary
    + gg.theme(axis_text_y=gg.element_blank())
)

output_file = pathlib.Path(f"{output_fig_dir}/feature_select_summary.png")
feature_select_gg.save(output_file, dpi=dpi, height=4, width=6)

feature_select_gg

