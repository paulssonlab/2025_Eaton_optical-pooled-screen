## crispri Container Instructions

### Requirements
A Linux environment with **Apptainer** installed

To install Apptainer in a linux kernel, run the following commands:

```
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

### Directory layout
Place the container, notebook, and data in the same folder:

```
project/
  crispri.sif
  Example_Analysis_Notebook.ipynb
  Example_Analysis_Notebook_Interactive.ipynb
  Example_Raw_Data/
```

File sources:
- `crispri.sif`: container from Zenodo (`doi:10.5281/zenodo.14537796`) at `/Containers/crispri.sif`
- `Example_Raw_Data/`: imaging data from Zenodo (`doi:10.5281/zenodo.14537796`)
- `Example_Analysis_Notebook.ipynb` and `Example_Analysis_Notebook_Interactive.ipynb`: from the GitHub repo `paulssonlab/2025_Eaton_optical-pooled-screen` at `/Image_Analysis/Example_Analysis`

### Run JupyterLab from the container
From inside `project/` (the directory containing `crispri.sif`), run:

```bash
apptainer exec --bind "$PWD" crispri.sif python -m jupyterlab --no-browser
```

This bind-mounts your current directory into the container, so the notebook and data remain accessible at the same relative paths.

You may follow the instructions in `Example_Analysis_Notebook.ipynb` to run an example analysis using pre-set parameters or use `Example_Analysis_Notebook_Interactive.ipynb` to run an analysis with custom parameters.

### Notes
- The data you want to analyze should be in a subdirectory of the folder that contains `crispri.sif` (or otherwise within the directory you bind).