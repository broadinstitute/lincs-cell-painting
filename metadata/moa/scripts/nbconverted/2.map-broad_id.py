#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# Samples and Drugs files from each version of the repurposing hub are read

# In[2]:


samples_2017 = pd.read_csv('clue/repurposing_samples_20170327.txt',usecols=['broad_id','pert_iname','InChIKey'],delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2017 = pd.read_csv('clue/repurposing_drugs_20170327.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')
samples_2018 = pd.read_csv('clue/repurposing_samples_20180907.txt', usecols=['broad_id','pert_iname','InChIKey'], delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2018 = pd.read_csv('clue/repurposing_drugs_20180907.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')
samples_2020 = pd.read_csv('clue/repurposing_samples_20200324.txt', usecols=['broad_id','pert_iname','InChIKey'], delimiter='\t', comment='!', encoding='iso-8859-1')
drugs_2020 = pd.read_csv('clue/repurposing_drugs_20200324.txt', usecols=['pert_iname','moa','target'], delimiter='\t', comment='!', encoding='iso-8859-1')


# In[3]:


# Maps the samples to their moa and target annotations

def target_annotation(samples, drugs):
    samples = samples.merge(drugs, on='pert_iname', how='left')
    return samples

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


# In[4]:


samples_2017 = target_annotation(samples_2017, drugs_2017)
samples_2018 = target_annotation(samples_2018, drugs_2018)
samples_2020 = target_annotation(samples_2020, drugs_2020)

samples_2017.head()


# In[5]:


samples_2017 = id_cleanup(samples_2017, "2017")
samples_2018 = id_cleanup(samples_2018, "2018")
samples_2020 = id_cleanup(samples_2020, "2020")

samples_2017.head()


# In[6]:


samples_2017 = group_by_InChIKey14(samples_2017, "2017")
samples_2018 = group_by_InChIKey14(samples_2018, "2018")
samples_2020 = group_by_InChIKey14(samples_2020, "2020")

samples_2017.head()


# The three dataframes are merged on InChIKey14

# In[7]:


merged_df = samples_2017.merge(samples_2018, on='InChIKey14', how='outer')
merged_df = merged_df.merge(samples_2020, on='InChIKey14', how='outer')
merged_df= merged_df.replace(to_replace='', value=np.nan)

merged_df.head()


# In[8]:


merged_df.to_csv('clue/broad_id_map.csv', index=False)


# Check if deprecated_broad_id has additional information not captured already

# In[9]:


deprecated_broad_id_2018 = pd.read_csv('clue/repurposing_samples_20180907.txt',usecols=['deprecated_broad_id','InChIKey'],delimiter='\t', comment='!', encoding='iso-8859-1')
deprecated_broad_id_2020 = pd.read_csv('clue/repurposing_samples_20200324.txt',usecols=['deprecated_broad_id','InChIKey'],delimiter='\t', comment='!', encoding='iso-8859-1')


# In[10]:


# Clean up depricated_broad_id

def dep_id_cleanup(df, year):
    df.deprecated_broad_id = df.deprecated_broad_id.apply(lambda x: str(x)[:13])
    df.InChIKey = df.InChIKey.apply(lambda x: str(x)[:14])
    df = df.drop_duplicates(['InChIKey','deprecated_broad_id']).reset_index(drop=True)
    df = df.rename(columns={'InChIKey':'InChIKey14',
                            'deprecated_broad_id':'deprecated_broad_id_'+year})

    return df

def dep_group_by_InChIKey14(df, year):
    df = df.groupby('InChIKey14')['deprecated_broad_id_'+year].apply(lambda x: '|'.join(np.unique(x))).reset_index()

    return df


# In[11]:


deprecated_broad_id_2018.dropna(inplace=True)
deprecated_broad_id_2018.reset_index(drop=True, inplace=True)

deprecated_broad_id_2020.dropna(inplace=True)
deprecated_broad_id_2020.reset_index(drop=True, inplace=True)


# In[12]:


deprecated_broad_id_2018 = dep_id_cleanup(deprecated_broad_id_2018, "2018")
deprecated_broad_id_2020 = dep_id_cleanup(deprecated_broad_id_2020, "2020")


# In[13]:


deprecated_broad_id_2018 = dep_group_by_InChIKey14(deprecated_broad_id_2018, "2018")
deprecated_broad_id_2020 = dep_group_by_InChIKey14(deprecated_broad_id_2020, "2020")


# In[14]:


unique_to_2018 = deprecated_broad_id_2018.merge(samples_2017, on='InChIKey14', how='left', indicator=True)._merge!='both'
unique_to_2020 = deprecated_broad_id_2020.merge(samples_2018, on='InChIKey14', how='left', indicator=True)._merge!='both'


# The following compounds are not present in the 2017 version

# In[15]:


print(deprecated_broad_id_2018[unique_to_2018].reset_index(drop=True).to_markdown())


# The following compounds are not present in the 2018 version

# In[16]:


print(deprecated_broad_id_2020[unique_to_2020].reset_index(drop=True).to_markdown())

