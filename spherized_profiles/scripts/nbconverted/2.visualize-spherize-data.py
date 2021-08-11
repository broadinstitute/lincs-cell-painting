#!/usr/bin/env python
# coding: utf-8

# ## Visualize a representation of the spherized LINCS Cell Painting dataset

# In[1]:


import umap
import pathlib
import numpy as np
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


np.random.seed(9876)


# In[3]:


profile_path = pathlib.Path("profiles")

batches = ["2016_04_01_a549_48hr_batch1", "2017_12_05_Batch2"]
norm_methods = ["whole_plate", "dmso"]

file_filler = "_dmso_spherized_profiles_with_input_normalized_by_"

output_dir = pathlib.Path("figures")
output_dir = {batch: pathlib.Path(output_dir, batch) for batch in batches}


# In[4]:


# Identify UMAP embeddings for all spherized profiles
embeddings = {batch: {} for batch in batches}
for batch in batches:
    for norm_method in norm_methods:
        file = pathlib.Path(profile_path, f"{batch}{file_filler}{norm_method}.csv.gz")
        print(f"Now obtaining UMAP embeddings for {file}...")
        # Load spherized data
        spherized_df = pd.read_csv(file)
        
        # Extract features
        cp_features = infer_cp_features(spherized_df)
        meta_features = infer_cp_features(spherized_df, metadata=True)

        # Fit UMAP
        reducer = umap.UMAP(random_state=123)
        embedding_df = reducer.fit_transform(spherized_df.loc[:, cp_features])

        embedding_df = pd.DataFrame(embedding_df)
        embedding_df.columns = ["UMAP_0", "UMAP_1"]
        embedding_df = pd.concat(
            [
                spherized_df.loc[:, meta_features],
                embedding_df
            ],
            axis="columns"
        )
        embedding_df = embedding_df.assign(dmso_label="DMSO")
        embedding_df.loc[embedding_df.Metadata_broad_sample != "DMSO", "dmso_label"] = "compound"

        embeddings[batch][norm_method] = embedding_df
        print("done.\n")


# ### Output a series of visualizations for the spherized profiles
# 
# 1. Batch 1 - Both normalization methods - UMAP highlighting compound vs. non-compound
# 2. Batch 1 - Both normalization methods - UMAP highlighting plate
# 3. Batch 2 - Both normalization methods - UMAP highlighting compound vs. non-compound
# 4. Batch 2 - Both normalization methods - UMAP highlighting plate
# 5. Batch 2 - Both normalization methods - UMAP highlighting different cell lines
# 6. Batch 2 - Both normalization methods - UMAP highlighting different time points
# 
# There will be a total of 12 figures, distributed in two pdf files.

# In[5]:


batch = "2016_04_01_a549_48hr_batch1"

plotlist = []
for norm_method in norm_methods:
    for color_type in ["Metadata_broad_sample", "Metadata_Plate"]:
        output_file = pathlib.Path(output_dir[batch], f"{batch}_{norm_method}_colorby{color_type}.png")
        output_file.parent.mkdir(exist_ok=True)
        
        label = f"Batch 1: Normalized by {norm_method.upper()}\nColored by {color_type}"

        embedding_gg = (
            gg.ggplot(embeddings[batch][norm_method], gg.aes(x="UMAP_0", y="UMAP_1"))
            + gg.geom_point(gg.aes(color=color_type), size=0.1, alpha=0.2)
            + gg.facet_grid("~dmso_label")
            + gg.ggtitle(label)
            + gg.theme_bw()
            + gg.xlab("UMAP X")
            + gg.ylab("UMAP Y")
            + gg.theme(
                legend_position="none",
                strip_text=gg.element_text(size=5),
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
                axis_text=gg.element_text(size=6),
                axis_title=gg.element_text(size=7),
                title=gg.element_text(size=7),
                figure_size=(5.5, 3)
            )
        )
        plotlist.append(embedding_gg)

output_file = pathlib.Path(output_dir[batch], f"{batch}_UMAPs.pdf")
output_file.parent.mkdir(exist_ok=True)
gg.save_as_pdf_pages(plotlist, output_file)


# In[6]:


batch = "2017_12_05_Batch2"

plotlist = []
for norm_method in norm_methods:
    for color_type in [
        "Metadata_broad_sample", "Metadata_Plate", "Metadata_cell_line", "Metadata_time_point"
    ]:
        output_file = pathlib.Path(output_dir[batch], f"{batch}_{norm_method}_colorby{color_type}.png")
        output_file.parent.mkdir(exist_ok=True)
        
        label = f"Batch 2: Normalized by {norm_method.upper()}\nColored by {color_type}"

        embedding_gg = (
            gg.ggplot(embeddings[batch][norm_method], gg.aes(x="UMAP_0", y="UMAP_1"))
            + gg.geom_point(gg.aes(color=color_type), size=0.1, alpha=0.2)
            + gg.facet_grid("~dmso_label")
            + gg.ggtitle(label)
            + gg.theme_bw()
            + gg.xlab("UMAP X")
            + gg.ylab("UMAP Y")
            + gg.theme(
                legend_position="none",
                strip_text=gg.element_text(size=5),
                strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
                axis_text=gg.element_text(size=6),
                axis_title=gg.element_text(size=7),
                title=gg.element_text(size=7),
                figure_size=(5.5, 3)
            )
        )
        
        if color_type in ["Metadata_cell_line", "Metadata_time_point"]:
            embedding_gg = (
                embedding_gg
                + gg.theme(
                    legend_position="right",
                    figure_size=(6, 3)
                )
            )
        
        plotlist.append(embedding_gg)

output_file = pathlib.Path(output_dir[batch], f"{batch}_UMAPs.pdf")
output_file.parent.mkdir(exist_ok=True)
gg.save_as_pdf_pages(plotlist, output_file)

