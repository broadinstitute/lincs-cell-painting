#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import rdkit.Chem.inchi as inchi


# Samples and Drugs files from each version of the repurposing hub are read

# In[2]:


repurposing_samples_2017 = pd.read_csv('clue/repurposing_samples_20170327.txt',usecols=['broad_id','pert_iname','InChIKey'],delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2017 = pd.read_csv('clue/repurposing_drugs_20170327.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')
repurposing_samples_2018a = pd.read_csv('clue/repurposing_samples_20180516.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2018a = pd.read_csv('clue/repurposing_drugs_20180516.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')
repurposing_samples_2018b = pd.read_csv('clue/repurposing_samples_20180907.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2018b = pd.read_csv('clue/repurposing_drugs_20180907.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')
repurposing_samples_2020 = pd.read_csv('clue/repurposing_samples_20200324.txt', usecols=['broad_id','pert_iname','InChIKey','deprecated_broad_id'], delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2020 = pd.read_csv('clue/repurposing_drugs_20200324.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')


# In[3]:


repurposing_samples_2017.dropna(subset=['InChIKey'], inplace=True)
repurposing_samples_2018a.dropna(subset=['InChIKey'], inplace=True)
repurposing_samples_2018b.dropna(subset=['InChIKey'], inplace=True)
repurposing_samples_2020.dropna(subset=['InChIKey'], inplace=True)

# 2017 version is missing deprecated_broad_id
repurposing_samples_2017['deprecated_broad_id'] = np.nan


# In[4]:


samples_2017 = repurposing_samples_2017.loc[:,['broad_id','pert_iname','InChIKey']]
deprecated_2017 = repurposing_samples_2017.loc[:,['deprecated_broad_id','pert_iname','InChIKey']]
samples_2018a = repurposing_samples_2018a.loc[:,['broad_id','pert_iname','InChIKey']]
deprecated_2018a = repurposing_samples_2018a.loc[:,['deprecated_broad_id','pert_iname','InChIKey']]
samples_2018b = repurposing_samples_2018b.loc[:,['broad_id','pert_iname','InChIKey']]
deprecated_2018b = repurposing_samples_2018b.loc[:,['deprecated_broad_id','pert_iname','InChIKey']]
samples_2020 = repurposing_samples_2020.loc[:,['broad_id','pert_iname','InChIKey']]
deprecated_2020 = repurposing_samples_2020.loc[:,['deprecated_broad_id','pert_iname','InChIKey']]


# In[5]:


# Merge broad_id and deprecated_broad_id

def merge_broad_dep_id(samples,deprecated):
    deprecated.dropna(subset=['deprecated_broad_id'], inplace=True)
    deprecated.rename(columns={'deprecated_broad_id':'broad_id'}, inplace=True)
    samples = samples.append(deprecated, sort=False, ignore_index=True)
    return samples

# Maps the samples to their moa and target annotations

def target_annotation(samples, drugs):
    samples = samples.merge(drugs, on='pert_iname', how='left')
    return samples


# In[6]:


# Replace InChI with InChIKey

def inchi_to_inchikey(df):
    df.InChIKey = df.InChIKey.apply(lambda x: inchi.InchiToInchiKey(x) if (x.startswith('InChI')) else x)
    return df

# The first 13 characters of the Broad ID and the first 14 characters of InChIKey are extracted
# Year names are appended to column names

def id_cleanup(df, year):
    df.broad_id = df.broad_id.apply(lambda x: str(x)[:13])
    df.InChIKey = df.InChIKey.apply(lambda x: str(x)[:14])
    df = df.drop_duplicates(['InChIKey','pert_iname','broad_id']).reset_index(drop=True)
    df = df.rename(columns={'pert_iname':'pert_iname_'+year,
                            'broad_id':'broad_id_'+year,
                            'InChIKey':'InChIKey14',
                            'moa':'moa_'+year,
                            'target':'target_'+year})
    return df

# Grouping samples using InChIKey14 while all other fields are pipe delimited

def group_by_InChIKey14(df, year):
    df = df.fillna('')
    df = df.groupby('InChIKey14').agg({'broad_id_'+year : lambda x: '|'.join(np.unique(x)),
                                       'pert_iname_'+year: lambda x: '|'.join(np.unique(x)),
                                       'moa_'+year : lambda x: merge_target(list(x)),
                                       'target_'+year: lambda x: merge_target(list(x))}).reset_index()
    return df

# This function deplicates target annotations and the final list is pipe delimited

def merge_target(target):
    out_target = '|'.join(np.unique(('|'.join(target)).split('|')))
    return out_target


# In[7]:


samples_2017 = merge_broad_dep_id(samples_2017, deprecated_2017)
samples_2018a = merge_broad_dep_id(samples_2018a, deprecated_2018a)
samples_2018b = merge_broad_dep_id(samples_2018b, deprecated_2018b)
samples_2020 = merge_broad_dep_id(samples_2020, deprecated_2020)


# In[8]:


samples_2017 = inchi_to_inchikey(samples_2017)
samples_2018a = inchi_to_inchikey(samples_2018a)
samples_2018b = inchi_to_inchikey(samples_2018b)
samples_2020 = inchi_to_inchikey(samples_2020)


# In[9]:


samples_2017 = target_annotation(samples_2017, drugs_2017)
samples_2018a = target_annotation(samples_2018a, drugs_2018a)
samples_2018b = target_annotation(samples_2018b, drugs_2018b)
samples_2020 = target_annotation(samples_2020, drugs_2020)

samples_2017.head()


# In[10]:


samples_2017 = id_cleanup(samples_2017, "2017")
samples_2018a = id_cleanup(samples_2018a, "2018a")
samples_2018b = id_cleanup(samples_2018b, "2018b")
samples_2020 = id_cleanup(samples_2020, "2020")

samples_2017.head()


# In[11]:


samples_2017 = group_by_InChIKey14(samples_2017, "2017")
samples_2018a = group_by_InChIKey14(samples_2018a, "2018a")
samples_2018b = group_by_InChIKey14(samples_2018b, "2018b")
samples_2020 = group_by_InChIKey14(samples_2020, "2020")

samples_2017.head()


# The four dataframes are merged on InChIKey14

# In[12]:


merged_df = samples_2017.merge(samples_2018a, on='InChIKey14', how='outer')
merged_df = merged_df.merge(samples_2018b, on='InChIKey14', how='outer')
merged_df = merged_df.merge(samples_2020, on='InChIKey14', how='outer')
merged_df= merged_df.replace(to_replace='', value=np.nan)

merged_df.head()


# In[13]:


merged_df.to_csv('clue/broad_id_map.csv', index=False)


# In[14]:


print('Total number of rows %d' % len(merged_df))


# The following compounds in the 2020 version have not been formatted correctly
# as InChI is being extracted instead of InChIKey for these compounds.

# In[15]:


print(merged_df.loc[merged_df.InChIKey14.str.startswith('InChI')][['InChIKey14','pert_iname_2020']].to_markdown())


