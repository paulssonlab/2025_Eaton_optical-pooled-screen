# Image_Analysis

This folder contains notebooks for mother machine time-lapse image analysis across multiple experimental datasets (organized by library/run ID). See the associated paper for detailed description of each library.

## Contents

#### GFP Library

- `lDE15/` - **two-variant GFP library (GFPmut2 vs DarkGFP)** used to benchmark end-to-end genotyping accuracy (FISH barcode calling + sequencing lookup).

#### CRISPRi Libraries

- `lDE20/` — **Essentialome-wide mismatch-CRISPRi screen**  pooled sgRNA library imaged in mother machine to quantify growth/morphology dynamics before and after CRISPRi induction.
- `lDE26/` — **HU-mCherry nucleoid-reporter CRISPRi library** selected strong sgRNA library enabling nucleoid phenotyping alongside standard growth/morphology.
- `lDE28/` — Same HU-mCherry CRISPRi library in a **ΔrelA background** (SpoT intact) used to dissect RelA-independent regulation in translation/growth-law phenotypes.
- `lDE30/` — Same HU-mCherry CRISPRi library in a **ΔrelA ΔspoT background** (no (p)ppGpp) used to test which scalings/phenotypes require (p)ppGpp.

#### Strain Isolates

- `lDE20_Validation/` — **Individually cloned CRISPRi isolates** (derived from library hits) imaged as single strains to validate pooled phenotypes.
- `DE32_Control/` — **Long-term control imaging dataset** (no CRISPRi perturbation) used to confirm stability of measured phenotype distributions over extended imaging.
- `Translation_Depletion_Strains/` — **Targeted translation-perturbation isolates** (initiation/elongation/translocation-related perturbations) used for mechanistic follow-up of size–growth and ribosome/(p)ppGpp-linked phenotypes.
- `Example_Analysis/` — **Small, self-contained example dataset** intended to demonstrate the analysis workflow end-to-end on a compact subset of data.