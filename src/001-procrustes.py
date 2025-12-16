# %%
import os
import numpy as np
from numpy.linalg import svd
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path as pth
from pyslingshot import Slingshot
import random
import anndata as ad
from sklearn.decomposition import PCA


# %% [markdown]
# ### Set Working Directory

# %%
os.chdir(pth.home() / 'traj-alignment')

# %% [markdown]
# ### Load Data

# %%

# set colnames
colnames =['Alcohol','Malic_acid','Ash','Alcalinity_of_ash','Magnesium','Total_phenols','Flavanoids','Nonflavanoid_phenols','Proanthocyanins','Color_intensity','Hue','OD280/OD315_of_diluted_wines','Proline']

colnames = [s.lower() for s in colnames]

wine = pd.read_csv(open('raw-data/wine.data', 'r'), delimiter = ',', index_col=None, names = colnames)
#%%
# wine_variety = wine.index.values-1
wine_variety = np.array(['x', 'y', 'z'])[wine.index.values-1]
wine_variety

# %%
wine.head()
# %% [markdown]
# ## Split Data
# Make two datasets from same original data. 
# - Overlap between samples but not features
# - Overlap between samples and features
# - No overlap between samples but overlap between features
# - No overlap between samples or features

# %%
# fix index
wine.index = range(0,wine.shape[0],1)

# Datasets with no sample or feature overlap 
n = wine.shape[0]
p = wine.shape[1]
n_rw = n // 2  
n_col = p // 2  

random.seed(1030)
rw_inds_a = random.sample(range(n), k=n_rw)
rw_inds_b = list(set(range(n)) - set(rw_inds_a))

a_variety = wine_variety[rw_inds_a]
b_variety = wine_variety[rw_inds_b]

col_inds_a = random.sample(range(p), k = n_col)
col_inds_b = list(set(range(p)) - set(col_inds_a))


# %%
a = wine.iloc[rw_inds_a, col_inds_a]
b = wine.iloc[rw_inds_b, col_inds_b]

#%% [markdown]
# #### Make AnnData

#%%
aobj = ad.AnnData(a)
bobj = ad.AnnData(b)

aobj.obs['variety'] = a_variety
bobj.obs['variety'] = b_variety

# aobj.cluster_labels=aobj.obs
# bobj.cluster_labels=bobj.obs

# %%
x = np.unique(np.concatenate([a.columns.values, b.columns.values]))
print('diff =', len(wine.columns.difference(x)),"columns")
# print(wine.columns)

# %% [markdown]
# ### Scale data and calculate pcs

# %%
# 0 mean unit variance
standardize = lambda x:  x.sub(x.mean(axis=0), axis = 1).div(x.std(axis =0), axis = 1)

#PCs
def get_PCs(X):
    pcs = PCA( n_components = 4).fit_transform(X)
    pcs = pd.DataFrame(pcs)
    pcs.columns = [f"PC{i}" for i in range(pcs.shape[1])]
    return pcs

aobj.layers["scld"]= standardize(a)
bobj.layers["scld"]= standardize(b)

aPCS = get_PCs(aobj.layers["scld"])
aPCS.index = aobj.obs_names
aobj.obsm["PCA"]= np.asarray(aPCS)

bPCS = get_PCs(bobj.layers["scld"])
bPCS.index = bobj.obs_names
bobj.obsm["PCA"]= np.asarray(bPCS)


# %% [markdown]
# ##### Plot PCS

# %%
print(aobj)
print(bobj)

# %%
a = (
    pd.DataFrame(aobj.obsm['PCA'])
    .assign(source = 'a')
)
b = (
    pd.DataFrame(bobj.obsm['PCA'])
    .assign(source = 'b')
)

plot_df = pd.concat([a, b], axis=0, ignore_index=True)
plot_df["variety"] = np.concatenate([a_variety, b_variety])

plt.figure(figsize=(6, 5))

for s in plot_df["source"].unique():
    mask = plot_df["source"] == s
    plt.scatter(
        plot_df.loc[mask, 0],
        plot_df.loc[mask, 1],
        c=plot_df.loc[mask, "variety"].map({"x": 0, "y": 1, "z": 2}),
        alpha=0.7
    )
    plt.xlabel("PC0")
    plt.ylabel("PC1")
    plt.legend(title="Wine variety")
    plt.show()


# %% [markdown]
# ### Calculate Trajectories

#%% 
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
custom_xlim = (-12, 12)
custom_ylim = (-12, 12)
# plt.setp(axes, xlim=custom_xlim, ylim=custom_ylim)

slingshot = Slingshot(aobj, 
                      celltype_key="variety", 
                      obsm_key="PCA", 
                      start_node=0, 
                      is_debugging="verbose")

slingshot.fit(num_epochs=1, debug_axes=axes)
plt.show()

#%%
fig, axes = plt.subplots(ncols=2, figsize=(12, 4))
axes[0].set_title("Clusters")
axes[1].set_title("Pseudotime")
slingshot.plotter.curves(axes[0], slingshot.curves)
slingshot.plotter.clusters(axes[0], labels=np.arange(slingshot.num_clusters), s=4, alpha=0.5)
slingshot.plotter.clusters(axes[1], color_mode="pseudotime", s=5)

# %% [markdown]
# ### Compute svd

# %%
aU, aS, aV = svd(a_scld)
bU, bS, bV = svd(b_scld)

print("Row Vectors = ", aU.shape)
print("Column Vectors = ",aS.shape)
print("Singular Values = ", aV.shape)

# %% [markdown]
# ### Orthogonal Procrustes alignment (map  __B__ --> __A__)
#
# Given centered and scaled matrices
# $$
# A, B \in \mathbb{R}^{n \times p},
# \qquad
# R \in \mathbb{R}^{p \times p},
# $$
# we estimate an orthogonal transformation that aligns __B__ to __A__ by solving
#
# $$
# R^* \;=\;
# \arg\min_{R \in \mathbb{R}^{p \times p}}
# \;\| B R - A \|_F^2
# \quad \text{subject to} \quad
# R^\top R = I .
# $$
#
# Let the cross-covariance matrix be
#
# $$
# M = B^\top A .
# $$
#
# If the singular value decomposition of \(M\) is
#
# $$
# M = U \Sigma V^\top ,
# $$
#
# then the optimal rotation matrix is given by
#
# $$
# R^* = U V^\top .
# $$
#



# %%
