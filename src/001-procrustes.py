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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% 
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as plt
import Path as pth

# %% [markdown]
# ### Set Working Directory

# %%
os.chdir('~/Documents/personal/R-stuff/traj-alignment')
# %% [markdown]
# ### Load Data

# %% 
wine = pd.read_csv(open('raw-data/wine.data', 'r'), delimiter = ',')

# %%
