#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import rdkit.Chem.inchi as inchi
import os.path


# Sample files from each version of the repurposing hub are read

# In[2]:


samples_2017 = pd.read_csv('clue/repurposing_samples_20170327.txt',usecols=['broad_id','pert_iname','InChIKey'],delimiter='\t', comment='!', encoding='iso-8859-1')
samples_2018a = pd.read_csv('clue/repurposing_samples_20180516.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')
samples_2018b = pd.read_csv('clue/repurposing_samples_20180907.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')
samples_2020 = pd.read_csv('clue/repurposing_samples_20200324.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')


# In[3]:


# 2017 version is missing deprecated_broad_id
samples_2017 = samples_2017.assign(deprecated_broad_id=np.nan)


# In[4]:


# The first 13 characters of the Broad ID and the first 14 characters of InChIKey are extracted
# Year names are appended to column names
# For some rows InChI are extracted instead of InChIKey -> these are fixed

def id_cleanup(df, year):
    df.dropna(subset=['InChIKey'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.InChIKey = df.InChIKey.apply(lambda x: inchi.InchiToInchiKey(x) if (x.startswith('InChI')) else x)
    df.broad_id = df.broad_id.apply(lambda x: str(x)[:13])
    df.InChIKey = df.InChIKey.apply(lambda x: str(x)[:14])
    df = df.drop_duplicates(['InChIKey','pert_iname','broad_id','deprecated_broad_id']).reset_index(drop=True)
    df = df.rename(columns={'pert_iname':'pert_iname_'+year,
                            'broad_id':'broad_id_'+year,
                            'deprecated_broad_id':'deprecated_broad_id_'+year,
                            'InChIKey':'InChIKey14'})
    return df

# Grouping samples using InChIKey14 while all other fields are pipe delimited

def group_by_InChIKey14(df, year):
    df = df.fillna('')
    df = df.groupby('InChIKey14').agg({'broad_id_'+year : lambda x: '|'.join(np.unique(x)),
                                       'deprecated_broad_id_'+year: lambda x: create_pipe_separated_list(list(x)),
                                       'pert_iname_'+year: lambda x: '|'.join(np.unique(x))}).reset_index()
    return df

# This function cleans up empty depredicated_broad_ids, de-duplicates them and then separates them with pipes

def create_pipe_separated_list(target):
    joined_target = ('|'.join(target)).split('|')
    while '' in joined_target:
        joined_target.remove('')
    out_target = '|'.join(np.unique(joined_target))
    return out_target


# In[5]:


samples_2017 = id_cleanup(samples_2017, "2017")
samples_2018a = id_cleanup(samples_2018a, "2018a")
samples_2018b = id_cleanup(samples_2018b, "2018b")
samples_2020 = id_cleanup(samples_2020, "2020")

samples_2017.head()


# In[6]:


samples_2017 = group_by_InChIKey14(samples_2017, "2017")
samples_2018a = group_by_InChIKey14(samples_2018a, "2018a")
samples_2018b = group_by_InChIKey14(samples_2018b, "2018b")
samples_2020 = group_by_InChIKey14(samples_2020, "2020")

samples_2017.head()


# The four data frames are merged on InChIKey14

# In[7]:


merged_df = samples_2017.merge(samples_2018a, on='InChIKey14', how='outer')
merged_df = merged_df.merge(samples_2018b, on='InChIKey14', how='outer')
merged_df = merged_df.merge(samples_2020, on='InChIKey14', how='outer')

merged_df= merged_df.replace(to_replace='', value=np.nan)

merged_df.head()


# In[8]:


map_file = os.path.join('clue', 'broad_id_map.csv')
merged_df.to_csv(map_file, index=False)


# In[9]:


print(merged_df.shape)

