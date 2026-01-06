# 2025_Eaton_optical-pooled-screen

This repository contains analysis code, data processing pipelines, and figures for the MARLIN (Multiplexed Assignment of RNA‐barcoded LINeages) optical pooled screening platform presented in Eaton et al., 2025.

Note that most analyses require the accompanying datasets deposited on Zenodo (`doi:10.5281/zenodo.14537796`).

## Repository Structure

- `Library_Design/` — Library design resources (e.g., sgRNA/barcode design and related analysis/code).
- `Agar_Pad_Image_Analysis/` — Notebooks/scripts for agar pad analysis (segmentation, phenotype extraction, downstream summaries).
- `Image_Analysis/` — Mother machine + optical pooled imaging analysis workflows (segmentation/tracking, barcode/FISH decoding, dataset-specific notebooks).
- `Sequencing/` — Sequencing processing and analysis for MARLIN libraries (e.g., nanopore and related summary outputs).
- `Replication_Runout/` — Notebooks for replication runout measurements and downstream analysis (e.g., deep sequencing / flow cytometry–based analyses).
- `Figures/` — Figure-generation notebooks and scripts (Figures 1–7; Extended Data Figs. 1–13), typically run against processed outputs from the Zenodo deposit.
- `containers/` — Apptainer/Singularity definitions and usage notes for reproducible execution of key analysis environments.
- `repros/` — Reproducibility helpers (e.g., dependency manifests/lockfiles and environment exports used to recreate analysis environments).
- `scripts/` — Small utility script to generate contents of `repros/`.
- `README.md` — This file.

  
## Instructions for Reproducing Figures

Manuscript figures (**Figs. 1–7** and **Extended Data Figs. 1–13**) are generated primarily from the Jupyter notebooks in `Figures/`, using processed tables and intermediate outputs deposited on Zenodo.

1. **Get the processed data**
   - Download the Zenodo deposit and unzip it locally.
   - Most figure notebooks expect to read *processed* outputs (feature tables, decoded barcode calls, summary statistics, etc.) from this archive.

2. **Run the figure notebooks**
   - Open the relevant notebooks in `Figures/` (typically organized by figure / panel).
   - Update the data root/path variables at the top of each notebook to point to your unzipped Zenodo folder.
   - Execute the notebook to regenerate the plots/panels; outputs are written to the figure output locations defined in each notebook.

## Instructions for Reproducing Agar Pad Analysis

To reproduce the agar pad image analysis end-to-end, use the container/setup guide for environment + Jupyter launch, and the example notebook for the full processing workflow.

- **1) Set up the environment + launch Jupyter**
  - Follow `containers/agarpad_instructions.md` to install Apptainer, place the container + notebook + example data in a single working directory, and start JupyterLab from inside the container.

- **2) Run the example notebook**
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
  - Follow `containers/crispri_instructions.md` to install Apptainer, place the container + notebook + example data in a single working directory, and start JupyterLab from inside the container.

- **2) Run the example notebook**
  - Open `Image_Analysis/Example_Analysis/Example_Analysis_Notebook.ipynb` and execute top-to-bottom (or use the interactive variant at `Image_Analysis/Example_Analysis/Example_Analysis_Notebook_Interactive.ipynb` if you want to tune parameters).
  - The notebook walks through the standard mother machine pipeline, including file extraction/preprocessing, trench detection/cropping, lineage/kymograph generation, segmentation + tracking, downstream quantification, and exporting final analysis tables/outputs back into the experiment folder.

## Instructions for Reproducing Library Sequencing

Raw Nanopore and NGS sequencing reads are available via NCBI BioProject **PRJNA1205775**.

Nanopore processing for this study is provided as a **Snakemake** pipeline in `sequencing/`. The workflow takes raw FAST5 files through basecalling and barcode assignment to produce a final **barcode → consensus/reference** table (`output.tsv`), with an accompanying QC/export notebook.

At a high level, the workflow chunks FAST5 files, basecalls with Guppy, aligns reads to the barcode graph (GraphAligner), assigns/group reads by barcode, calls per-barcode consensus sequences (Medaka), and writes a merged `output.tsv` suitable for downstream mapping/QC.

For example, to reproduce the analysis for `lDE26`:
- See `Sequencing/README.md` for an overview, required inputs, and an example SLURM Snakemake command.
- Use `Sequencing/lDE26/lDE26_Sequencing.smk` + `Sequencing/lDE26/configs/` + `Sequencing/lDE26/reference_sequences/` to run the pipeline.
- Use `Sequencing/lDE26/Barcode_QC_and_Export.ipynb` for QC and exporting summary tables.