## agarpad Container Instructions

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
  agarpad.sif
  Example_Analysis_Notebook.ipynb
  Example_Analysis_Notebook_Interactive.ipynb
  Example_Raw_Data/
```

File sources:
- `agarpad.sif`: container from Zenodo (`doi:10.5281/zenodo.14537795`) at `/Containers/agarpad.sif`
- `Example_Raw_Data/`: imaging data from Zenodo (`doi:10.5281/zenodo.14537795`)
- `Example_Analysis_Notebook.ipynb` and `Example_Analysis_Notebook_Interactive.ipynb`: from the GitHub repo `paulssonlab/2025_Eaton_optical-pooled-screen` at `/Agar_Pad_Image_Analysis/Example_Analysis`

### Run JupyterLab from the container
From inside `project/` (the directory containing `agarpad.sif`), run:

```bash
apptainer exec --bind "$PWD" agarpad.sif python -m jupyterlab --no-browser
```

Or, if you would like to use a GPU to accelerate segmentation:

```bash
apptainer exec --bind "$PWD" --nv agarpad.sif python -m jupyterlab --no-browser
```

This bind-mounts your current directory into the container, so the notebook and data remain accessible at the same relative paths.

You may follow the instructions in `Example_Analysis_Notebook.ipynb` to run an example analysis using pre-set parameters or use `Example_Analysis_Notebook_Interactive.ipynb` to run an analysis with custom parameters.

### Notes
- The data you want to analyze should be in a subdirectory of the folder that contains `agarpad.sif` (or otherwise within the directory you bind).