#!/usr/bin/env python
# coding: utf-8

# In[26]:


import umap
import pathlib
import pandas as pd
import plotnine as gg

from pycytominer.cyto_utils import infer_cp_features
from cytominer_eval import evaluate


# In[4]:


batches = ["2016_04_01_a549_48hr_batch1", "2017_12_05_Batch2"]
suffixes = ["_spherized_profiles_dmso.csv.gz", "_spherized_profiles_whole_plate.csv.gz"]
profile_path = pathlib.Path("profiles")


# In[5]:


for batch in batches:
    for suffix in suffixes:
        spherized_file = pathlib.Path(profile_path, f"{batch}{suffix}")
        spherized_df = pd.read_csv(spherized_file, low_memory=False)


# In[ ]:


evaluate(
    spherized_df,
    features=cp_features,
    meta_features=meta_features,
    replicate_groups=["Metadata_cell_line", "Metadata_broad_sample", "Metadata_mg_per_ml"]
)


# In[28]:


spherized_df.head()


# In[ ]:





# In[ ]:





# In[21]:


cp_features = infer_cp_features(spherized_df)
meta_features = infer_cp_features(spherized_df, metadata=True)

# Apply UMAP
reducer = umap.UMAP(random_state=123)
embedding_df = reducer.fit_transform(spherized_df.loc[:, cp_features])

embedding_df = pd.DataFrame(embedding_df)
embedding_df.columns = ["umap_x", "umap_y"]
embedding_df = pd.concat(
    [
        spherized_df.loc[:, meta_features],
        embedding_df
    ],
    axis="columns"
)
embedding_df = embedding_df.assign(dmso_label="DMSO")
embedding_df.loc[embedding_df.Metadata_broad_sample != "DMSO", "dmso_label"] = "compound"


# In[22]:


embedding_df.head()


# In[24]:


embedding_df.dtypes


# In[23]:


(
    gg.ggplot(embedding_df, gg.aes(x="umap_x", y="umap_y"))
    + gg.geom_point()
)


# In[25]:


embedding_df.plot(x="umap_x", y="umap_y", kind="scatter")


# In[15]:


(
           gg.ggplot(embedding_df, gg.aes(x="x", y="y")) +
           gg.geom_point(gg.aes(size="Metadata_mg_per_ml", color="Metadata_broad_sample"), alpha=0.5) +
           gg.facet_grid("~dmso_label") +
           gg.theme_bw() +
           gg.theme(legend_position="none", strip_background=gg.element_rect(colour="black", fill="#fdfff4"))
       )


# In[ ]:




