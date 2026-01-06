# 2025_Eaton_optical-pooled-screen

This repository contains analysis code, data processing pipelines, and figures for the MARLIN (Multiplexed Assignment of RNA‐barcoded LINeages) optical pooled screening platform presented in Eaton et al., 2025.

## Repository Structure

- **Library_Design/**: Code and protocols for designing the mismatch‐CRISPRi library targeting 585 essential genes with ~29,738 sgRNAs and unique 30-bit FISH barcodes.
- **Agar_Pad_Image_Analysis/**: Jupyter notebooks and scripts for segmenting and extracting timelapse phenotypes (length, width, growth rate, etc) from agar pad images.
- **Image_Analysis/**: Custom pipelines for analyzing MARLIN imaging data including decoding combinatorial FISH barcodes and high throughput segmentation/tracking of mother machine data. 
- **Sequencing/**: Scripts for processing nanopore sequencing data for MARLIN libraries.
- **Replication_Runout/**: Analysis notebooks for processing deep sequencing and flow cytometry data used to measure replication defects.
- **Figures/**: Scripts and notebooks to generate manuscript figures.
- **README.md**: This file.

## Version Text TBD

## Instructions for Reproducing Figures

- Reference the Zonodo archive
- Find all instances where there is a file in the zenodo referenced and make the path ./zipfolder_name/path/to/file
- List all instances where raw data that was not deposited is used and say available upon request

## Instructions for Reproducing Agar Pad Analysis

To reproduce the agar pad image analysis end-to-end, use the container/setup guide for environment + Jupyter launch, and the example notebook for the full processing workflow. Note that you will also need access to the data deposited at Zonodo (`doi:10.5281/zenodo.14537796`).

- **1) Set up the analysis environment + launch Jupyter**
  - Follow `containers/agarpad_instructions.md` to install Apptainer, place the container + notebook + example data in a single working directory, and start JupyterLab from inside the container.

- **2) Run the example workflow notebook**
  - Open `Agar_Pad_Image_Analysis/Example_Analysis/Example_Analysis_Notebook.ipynb` and execute top-to-bottom.
  - You will (i) point the notebook at the example dataset (e.g., `./Example_Raw_Data/Agar_Pad`), and choose **local vs SLURM** execution (via `analysis_is_local` and `dask_workingdir`).
  - The notebook then performs:
    - **Flat-field generation** (dark + fluorescence channels)
    - **File conversion + extraction**: ND2 → TIFF → HDF5 (via TrenchRipper extractors)
    - **Segmentation**: phase segmentation with **Omnipose** (GPU recommended/required), plus **nucleoid segmentation** (cell-wise Otsu)
    - **Feature extraction + merge**: regionprops tables → a single per-cell analysis table
  - Final outputs are written back into the experiment folder, including the regionprops parquet folders (e.g., `cell_analysis/`, `nuc_analysis/`) and a combined per-cell table (e.g., `Final_Output_df.pkl` and `Final_Output_df.csv`).

## Instructions for Reproducing Mother Machine Analysis

To reproduce the mother-machine analysis end-to-end, use the container/setup guide to launch Jupyter in a reproducible environment, then run the example notebook to execute the full TrenchRipper workflow.

- **1) Set up the environment + launch Jupyter**
  - Follow the container instructions to install Apptainer, organize your working directory (container + notebooks + example data), and start JupyterLab from inside the container.

- **2) Run the example mother-machine workflow notebook**
  - Open `Image_Analysis/Example_Analysis/Example_Analysis_Notebook.ipynb` and execute top-to-bottom (or use the interactive variant at `Image_Analysis/Example_Analysis/Example_Analysis_Notebook_Interactive.ipynb` if you want to tune parameters).
  - The notebook walks through the standard mother-machine pipeline, including file extraction/preprocessing, trench detection/cropping, lineage/kymograph generation, segmentation + tracking, downstream quantification, and exporting final analysis tables/outputs back into the experiment folder.

## Instructions for Reproducing Library Sequencing

## To-do

- Also make sure to deposit toy data in the zenodo (I'll have to see how to push the version number of the repo)
- I added the agar pad data and _??_
- Still need to crop one FOV from a mother machine dataset and test