# CRABS summer school 2025 - scRNA-seq practicals

## Installation

```bash
conda create -n crabs-py311 python=3.11 --yes && conda activate crabs-py311
pip install -e ".[dev,jupyter]"
pre-commit install

python -m ipykernel install --user --name crabs-py311 --display-name "crabs-py311"
```
