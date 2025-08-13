# %% [markdown]
# # Single-cell RNA sequencing data analysis tutorial
#
# In this first session, we will preprocess and cell type a [dataset of peripheral blood mononuclear cells](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0) (PBMCs).
#
# Useful links:
# * [anndata](https://anndata.readthedocs.io/en/stable/)
# * [scanpy](https://scanpy.readthedocs.io/en/stable/)
# * [Single-cell best practices book](https://www.sc-best-practices.org/preamble.html)

# %% [markdown]
# ## Library imports

# %%
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import mplscience
import seaborn as sns

import scanpy as sc

# `DATA_DIR` is a pathlib Path pointing to this repo's data directory.
# You can specify a path with the syntax `DATA_DIR / path / to / file`.
# Note: This import command only works if you followed the first installation instructions!
from crabs import DATA_DIR

# %% [markdown]
# ## General settings

# %%
# Set the verbosity of Scanpy to 2
sc.settings.verbosity = 2

# %% [markdown]
# ## Constants

# %% [markdown]
# ## Function definitions

# %% [markdown]
# ## Data loading

# %%
# Load the PBMC data into the AnnData format
adata = sc.read_10x_h5(DATA_DIR / "pbmc" / "raw" / "pbmc_10k_v3_filtered_feature_bc_matrix.h5")
adata

# %% [markdown]
# ## Data exploration

# %%
# Follow the instruction issued when loading the data
# Makes all gene names unique
adata.var_names_make_unique()

# %%
# How many cells and genes does the data contain?
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# %%
# Which format is the data saved as? Why?
print(f"Count data type: {type(adata.X)}")

# scRNA-seq data is very sparse, i.e., has many zeros entries
print(f"Percentage of non-zero counts: {adata.X.getnnz() / adata.n_obs / adata.n_vars * 100:.3}%")

# %% [markdown]
# ## Quality control

# %%
# Remove cells with less than 100 transcripts
sc.pp.filter_cells(adata, min_counts=100)
adata

# %%
# Remove genes expressed in less than 10 cells
sc.pp.filter_genes(adata, min_cells=10)
adata

# %%
# Detect mitochondrial genes and add a boolean flag to the gene metadata
adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT"))

# %%
# Detect ribosomal genes and add a boolean flag to the gene metadata
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

# %%
# Compute cell-wise quality metrics
# - Total counts
# - The number of genes expressed each cell
# - Total number of counts for mitochondrial genes
# - Proportion of total counts which are mitochondrial
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=[20], log1p=True, inplace=True)

# %%
# Plot the number of transcripts present in each cell against the percentage of observed mitochondrial genes
with mplscience.style_context():
    sns.set_style(style="whitegrid")

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot(data=adata.obs, x="total_counts", y="pct_counts_mt", s=1, color="grey", ax=ax)
    ax.set(xlabel="Total counts per cell", ylabel="Mitochondrial count percentage")

# %%
# Plot the number of transcripts present in each cell against the number of genes
with mplscience.style_context():
    sns.set_style(style="whitegrid")

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot(data=adata.obs, x="total_counts", y="n_genes_by_counts", s=1, color="grey", ax=ax)
    ax.set(xlabel="Total counts per cell", ylabel="Number of genes per cell")
    plt.show()

# %%
# Remove likely dead/dying cells and data outliers
outlier_obs_mask = (adata.obs["total_counts"] > 30000).values | (adata.obs["pct_counts_mt"] > 45).values

adata = adata[~outlier_obs_mask, :].copy()
adata

# %%
# Replot above's plots
with mplscience.style_context():
    sns.set_style(style="whitegrid")

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot(data=adata.obs, x="total_counts", y="pct_counts_mt", s=1, color="grey", ax=ax)
    ax.set(xlabel="Total counts per cell", ylabel="Mitochondrial count percentage")

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.scatterplot(data=adata.obs, x="total_counts", y="n_genes_by_counts", s=1, color="grey", ax=ax)
    ax.set(xlabel="Total counts per cell", ylabel="Number of genes per cell")
    plt.show()

# %% [markdown]
# ## Data preprocessing

# %% [markdown]
# ### Doublet detection

# %%
# Compute doublet statistics with scrublet
sc.pp.scrublet(adata)

# %%
# Visualize the distribution of the doublet score
# Above which threshold is a cell considered a doublet?
with mplscience.style_context():
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.kdeplot(data=adata.obs, x="doublet_score", color="gray", fill=True, ax=ax)
    ax.vlines(
        x=adata.obs.loc[adata.obs["predicted_doublet"], "doublet_score"].min(),
        ymin=0,
        ymax=ax.get_ylim()[1],
        colors="black",
        linestyles="--",
    )

    ax.set(xlabel="Doublet score")

    plt.show()

# %%
# How many cells does scrublet identify as putative doublets?
print(f"Number of cells identified as doublets: {adata.obs['predicted_doublet'].sum()}")
print(f"Percentage of cells identified as doublets: {adata.obs['predicted_doublet'].sum() / adata.n_obs * 100.:3}%")

# %% [markdown]
# ### Transformation and feature selection

# %%
# Normalize the cells to their median library size (total number of transcripts per cell)
sc.pp.normalize_total(adata)

# %%
# Log1p transform the data
# Why do we use log1p transformation?
sc.pp.log1p(adata)

# %%
# Select the 4000 most highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=False)

# %% [markdown]
# ### Dimensionality reduction

# %%
# Compute the principal component embedding
sc.tl.pca(adata)

# %%
# What is a suitable number of principle components to use

n_pcs = len(adata.uns["pca"]["variance_ratio"])
with mplscience.style_context():
    sns.set_style("whitegrid")

    df = pd.DataFrame({"variance_ratio": adata.uns["pca"]["variance_ratio"], "pc": np.arange(start=1, stop=n_pcs + 1)})

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(data=df, x="pc", y="variance_ratio", color="black", ax=ax)
    sns.scatterplot(data=df, x="pc", y="variance_ratio", color="black", ax=ax)
    ax.set(xlabel="Principal component", ylabel="Variance ratio")

    plt.show()

# %%
# Compute a k-nearest neighbor graph with k=30
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)

# %%
# Compute and plot the UMAP embedding of the data
sc.tl.umap(adata)

with mplscience.style_context():
    sc.pl.umap(adata)
    plt.show()

# %% [markdown]
# ## Data clustering

# %%
# Compute a Leiden clustering of the data
sc.tl.leiden(adata)

# %%
# Display the UMAP embedding, colored by leiden cluster
with mplscience.style_context():
    sc.pl.umap(adata, color="leiden", legend_loc="on data")
    plt.show()

# %%
# Is there a cluster of cells likely comprised of doublets?
with mplscience.style_context():
    sc.pl.umap(adata, color="doublet_score", legend_loc="on data", title="Doublet score")
    plt.show()

with mplscience.style_context():
    sns.set_style("whitegrid")
    palette = dict(zip(adata.obs["leiden"].cat.categories, adata.uns["leiden_colors"]))
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(data=adata.obs, x="leiden", y="doublet_score", hue="leiden", palette=palette, ax=ax)
    ax.set(xlabel="Leiden cluster", ylabel="Doublet score")
    ax.legend().remove()
    plt.show()

# %% [markdown]
# ## Cell typing

# %%
# Given the following cell type markers, assign cell types to clusters
# B cells: MS4A1, CD79A, IGHM
# T cells: CD4, CD8A, CD3D, CD3E, TRAC, LTB, IL7R, CCL4
# Natural killer cells: NKG7, KLRD1, PRF1, CD74, FCGR3A and low T cell marker expression
# Monocytes: CD14, CD163, TYROBP, LYZ, CSF1R, CD68

# %%
sc.tl.rank_genes_groups(adata=adata, groupby="leiden")

# %%
marker_dictionary = {
    "B": ["MS4A1", "CD79A", "IGHM"],
    "T": ["CD4", "CD8A", "CD3D", "CD3E", "TRAC", "LTB", "IL7R", "CCL4"],
    "NK": ["NKG7", "KLRD1", "PRF1", "CD74", "FCGR3A"],
    "Monocytes": ["CD14", "CD163", "TYROBP", "LYZ", "CSF1R", "CD68"],
}

with mplscience.style_context():
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby="leiden",
        var_names=marker_dictionary,
        standard_scale="var",
    )
