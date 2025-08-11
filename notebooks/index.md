# CRABS summer school 2025 - scRNA-seq practicals

## Installation

There are two options to complete our practical session:

1. You can clone this GitHub repo and install the corresponding package:

```bash
# If you use mamba, simply replace `conda` by `mamba`
conda create -n crabs-py311 python=3.11 --yes && conda activate crabs-py311

git clone https://github.com/WeilerP/crabs2025.git
cd crabs2025
pip install -e ".[jupyter]"

python -m ipykernel install --user --name crabs-py311 --display-name "crabs-py311"
```

2. Install the required packages

```bash
# If you use mamba, simply replace `conda` by `mamba`
conda create -n crabs-py311 python=3.11 --yes && conda activate crabs-py311

pip install anndata, igraph, leidenalg, scanpy, scikit-image, jupyterlab, ipywidget

python -m ipykernel install --user --name crabs-py311 --display-name "crabs-py311"
```

## Jupter lab

To run Jupyter Notebooks you can either use an IDE like [Visual Studio Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks), or [start a Jupyter server from your terminal](https://jupyterlab.readthedocs.io/en/stable/getting_started/starting.html).
