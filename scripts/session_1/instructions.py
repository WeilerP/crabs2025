# %% [markdown]
# # Single-cell RNA sequencing data analysis tutorial
#
# In this first session, we will preprocess and cell type a [dataset of peripheral blood mononuclear cells](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0) (PBMCs).
#
# If you cloned the GitHub repo and installed the corresponding package, save the data in the `data/` directory.
#
# Useful links:
# * [anndata](https://anndata.readthedocs.io/en/stable/)
# * [scanpy](https://scanpy.readthedocs.io/en/stable/)
# * [Single-cell best practices book](https://www.sc-best-practices.org/preamble.html)

# %% [markdown]
# ## Library imports

# %%
# `DATA_DIR` is a pathlib Path pointing to this repo's data directory.
# You can specify a path with the syntax `DATA_DIR / path / to / file`.
# Note: This import command only works if you followed the first installation instructions!
from crabs import DATA_DIR  # noqa

# %% [markdown]
# ## General settings

# %%
# Set the verbosity of Scanpy to 2

# %% [markdown]
# ## Constants

# %% [markdown]
# ## Function definitions

# %% [markdown]
# ## Data loading

# %%
# Load the PBMC data into the AnnData format


# %% [markdown]
# ## Data exploration

# %%
# Follow the instruction issued when loading the data


# %%
# How many cells and genes does the data contain?


# %%
# Which format is the data saved as? Why?

# scRNA-seq data is very sparse, i.e., has many zeros entries


# %% [markdown]
# ## Quality control

# %%
# Remove cells with less than 100 transcripts


# %%
# Remove genes expressed in less than 10 cells


# %%
# Detect mitochondrial genes and add a boolean flag to the gene metadata


# %%
# Detect ribosomal genes and add a boolean flag to the gene metadata


# %%
# Compute cell-wise quality metrics
# - Total counts
# - The number of genes expressed each cell
# - Total number of counts for mitochondrial genes
# - Proportion of total counts which are mitochondrial


# %%
# Plot the number of transcripts present in each cell against the percentage of observed mitochondrial genes


# %%
# Plot the number of transcripts present in each cell against the number of genes


# %%
# Remove likely dead/dying cells and data outliers


# %%
# Replot above's plots


# %% [markdown]
# ## Data preprocessing

# %% [markdown]
# ### Doublet detection

# %%
# Compute doublet statistics with scrublet


# %%
# Visualize the distribution of the doublet score
# Above which threshold is a cell considered a doublet?


# %%
# How many cells does scrublet identify as putative doublets?


# %% [markdown]
# ### Transformation and feature selection

# %%
# Normalize the cells to their median library size (total number of transcripts per cell)


# %%
# Log1p transform the data
# Why do we use log1p transformation?


# %%
# Select the 4000 most highly variable genes


# %% [markdown]
# ### Dimensionality reduction

# %%
# Compute the principal component embedding


# %%
# What is a suitable number of principle components to use


# %%
# Compute a k-nearest neighbor graph with k=30


# %%
# Compute and plot the UMAP embedding of the data


# %% [markdown]
# ## Data clustering

# %%
# Compute a Leiden clustering of the data


# %%
# Display the UMAP embedding, colored by leiden cluster


# %%
# Is there a cluster of cells likely comprised of doublets?


# %% [markdown]
# ## Cell typing

# %%
# Given the following cell type markers, assign cell types to clusters
# B cells: MS4A1, CD79A, IGHM
# T cells: CD4, CD8A, CD3D, CD3E, TRAC, LTB, IL7R, CCL4
# Natural killer cells: NKG7, KLRD1, PRF1, CD74, FCGR3A and low T cell marker expression
# Monocytes: CD14, CD163, TYROBP, LYZ, CSF1R, CD68

# %%

# %%
