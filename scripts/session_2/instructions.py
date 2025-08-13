# %% [markdown]
# # Single-cell RNA sequencing data analysis tutorial
#
# In this session, you will process another dataset of dataset of peripheral blood mononuclear cells (PBMCs), ultimately run into an issue, and get the chance to brainstorm about possible solutions to resolve the problem.
#
# You can download the needed data from [here](https://github.com/JinmiaoChenLab/Batch-effect-removal-benchmarking/tree/master/Data/dataset5); [this]() is the corresponding paper.
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
# Load the data of the first dataset (b1) into the AnnData format


# %%
# Load the data of the second dataset (b2) into the AnnData format

# %% [markdown]
# ## Data exploration

# %%
# Convert the count data into an appropriate, data-efficient format


# %%
# Does the cell metadata include information to uniquely identify each dataset?
# If not, add the dataset-specific identifiers `"sample_1"` and `"sample_2"` as
# the column `"sample"` to each AnnData's cell metadata (.obs).

# %%
# Combine the two AnnData objetcs into a single one, retaining all genes

# %%
# Delete the individual data objects to reduce the memory footprint

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
# ### Transformation and feature selection

# %%
# Normalize the cells to their median library size (total number of transcripts per cell)


# %%
# Log1p transform the data


# %%
# Select highly variable genes


# %% [markdown]
# ### Dimensionality reduction

# %%
# Compute the principal component embedding


# %%
# What is a suitbale number of principle components to use


# %%
# Compute a k-nearest neighbor graph with k=30


# %%
# Compute and plot the UMAP embedding of the data, colored by cell type and by sample
