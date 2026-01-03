## crispri Container Instructions

### Requirements
- A Linux environment with **Apptainer** installed (WSL2 is fine if Apptainer is installed inside the distro)

### Directory layout
Place the container, notebook, and data in the same folder:

```
project/
  crispri.sif
  Example_Analysis_Notebook.ipynb
  Example_Raw_Data/
```

File sources:
- `crispri.sif`: container from Zenodo (`doi:10.5281/zenodo.14537796`) at `/Containers/crispri.sif`
- `Example_Raw_Data/`: imaging data from Zenodo (`doi:10.5281/zenodo.14537796`)
- `Example_Analysis_Notebook.ipynb`: from the GitHub repo `paulssonlab/2025_Eaton_optical-pooled-screen`

### Run JupyterLab from the container
From inside `project/` (the directory containing `crispri.sif`), run:

```bash
apptainer exec --bind "$PWD" crispri.sif python -m jupyterlab --no-browser
```

This bind-mounts your current directory into the container, so the notebook and data remain accessible at the same relative paths.

### Notes
- The data you want to analyze should be in a subdirectory of the folder that contains `crispri.sif` (or otherwise within the directory you bind).
- If you run out of memory in WSL2, WSL may be limiting available RAM; increase the WSL memory limit if needed.
