# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: all,-execution,-papermill,-trusted
#     notebook_metadata_filter: -jupytext.text_representation.jupytext_version
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: traj-alignment
#     language: python
#     name: python3
# ---

# %%
import os
from pathlib import Path as pth
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from numpy.linalg import svd
from sklearn.decomposition import PCA
from pyslingshot import Slingshot
import anndata as ad
import scanpy as sc


# %% [markdown]
# ### Set Working Directory

# %%
os.chdir(pth.home() / 'dev/trajectory-alignment')

# %% [markdown]
# ### Load and Filter Data
# [link to data](https://covid19.cog.sanger.ac.uk/submissions/release2/meyer_nikolic_covid_pbmc_raw.h5ad)

# %%

file_path = pth('processed-data/covid-pbmc-filtered.h5ad')
if not file_path.exists() :
    adata = ad.read_h5ad('raw-data/meyer_nikolic_covid_pbmc_raw.h5ad')

    # filter for T CD4+
    adata = adata[adata.obs.annotation_broad == 'T CD4+']

    # filter for adults
    adata = adata[adata.obs.Age_group == 'Adult']

    # filter for Europeans
    adata = adata[adata.obs.Ethnicity == 'EUR']

    adata.write_h5ad(file_path)
    
else:
    adata = ad.read_h5ad(file_path)


# %%
import re

def snake_case(name: str) -> str:
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    s2 = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1)
    return s2.lower()

# rename obs 
adata.obs.columns = [snake_case(name) for name in adata.obs.columns]

adata.obs

# %% [markdown]
# #### Normalize 

# %%
# scanpy normalize_total and log1p functions
scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
# log1p transform - log the data and adds a pseudo-count of 1
scales_counts = sc.pp.log1p(scales_counts["X"], copy=True)

adata.layers['logcounts'] = scales_counts

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(adata.obs["n_count_rna"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(scales_counts.sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Shifted logarithm")
plt.show()

# %%
adata

# %% [markdown]
# #### Select Highest Variable Genes

# %%
# compute variable genes on logcounts
sc.pp.highly_variable_genes(adata, layer = 'logcounts', min_mean=0.0125, max_mean=3, min_disp=0.6)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))

#plot variable genes
sc.pl.highly_variable_genes(adata)


# %%
# subset for variable genes in the dataset
adata = adata[:, adata.var['highly_variable']]

# %% [markdown]
# #### Scale

# %%
# scales logcounts, writes to adata.X
sc.pp.scale(adata, layer= 'logcounts', max_value=10, zero_center = False )

# %% [markdown]
# #### Reduced Dimensions

# %%
adata

# %%
# pca on X written to X
sc.tl.pca(adata, svd_solver='arpack', n_comps=50, random_state = 129)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40)
sc.tl.umap(adata, random_state=42)


# %%
adata

# %%

sc.pl.embedding(adata, basis = 'X_umap', color='covid_severity')
sc.pl.embedding(adata, basis = 'X_umap', color='smoker')
sc.pl.embedding(adata, basis = 'X_pca', color='covid_severity')
sc.pl.embedding(adata, basis = 'X_pca', color='smoker')

# %%
