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