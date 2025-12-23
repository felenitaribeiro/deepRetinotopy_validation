# Repo with code for the validation of deepRetinotopy toolbox

This repository contains the code to validate the [deepRetinotopy toolbox](https://github.com/felenitaribeiro/deepRetinotopy_TheToolbox), accompanying the following publication:

[Add preprint here](https://www.biorxiv.org/content/10.1101/2025.11.27.690210v2)

## How to run the validation code

### 1. Create a conda environment

```bash
conda create -n deepretinotopy_validation python=3.8 ipykernel -y
conda activate deepretinotopy_validation
pip install -r requirements.txt;
```

### 2. Clone repos

```bash
git clone https://github.com/felenitaribeiro/deepRetinotopy_validation.git
cd deepRetinotopy_validation
git clone https://github.com/felenitaribeiro/deepRetinotopy_TheToolbox.git
```

### 3. Manuscript's analyses

You can find all code and plotting functionalities for reproducing the analyses reported in our manuscript in the *notebooks* folder, where notebooks are named according to the manuscript sections. We will continue to update the notebooks to enable reproducibility, as there are many hard-coded paths. However, one can already inspect the figures rendered in the notebooks to compare with our manuscript.