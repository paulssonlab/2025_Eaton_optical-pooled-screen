# Sequencing

This directory contains a Snakemake-based pipeline for processing Oxford Nanopore amplicon sequencing from MARLIN libraries, with the goal of producing a **barcode → variant/amplicon consensus** table (plus QC) that can be used downstream for optical pooled screening analyses.

The primary example workflow included here is **`lDE26`** (see `lDE26_Sequencing.smk` and the accompanying QC notebook).

If you trying to implement your own analysis, and not reproducing the sequencing results from the paper, consider instead using the better maintained repo at `https://github.com/paulssonlab/marlin-nanopore`

---

## Deposited sequencing data

Raw sequencing reads for Nanopore and NGS libraries have been deposited in the NCBI BioProject database under accession number **PRJNA1205775**.

## Directory layout

> (Exact filenames may vary by run, but the roles of each subdirectory are consistent.)

- `configs/`  
  Run-specific YAML configuration files (e.g., paths to input FAST5, scratch/output locations, flowcell/kit settings, chunk sizes, read-length filters, barcode depth thresholds, and subsampling settings).

- `reference_sequences/`  
  Reference sequence assets for alignment. For `lDE26`, the key input is a **GFA graph** used by GraphAligner.

- `scripts/`  
  Small Python utilities called by Snakemake rules/checkpoints.

- `profiles/`  
  Snakemake execution profiles (e.g., Slurm), used to run the workflow on an HPC cluster and request appropriate resources (notably a GPU for Guppy basecalling).

- `QC_outputs/`  
  Output folder used for QC plots/tables and exported summary data products.

---

## Example analysis: `lDE26`

### Running `lDE26`

1. Edit the config file in `configs/` (paths + run parameters).
2. Run Snakemake with your preferred execution mode (local or cluster profile). For cluster usage, use the provided Snakemake profile under `profiles/`.

The typical command to run this pipeline on `SLURM` is:

`snakemake -s lDE26_Sequencing.smk --profile profiles/slurm aggregate_outputs scratch_path/output.tsv`

After the workflow completes, the main product is the merged `output.tsv` written under the configured `scratch_path`.

### Required starting inputs

For each pipeline, you must define in the config file:

- **Raw Nanopore input**
  - `fast5_path`: A directory containing the input **`.fast5`** files.

- **Reference graph**
  - `ref_gfa_path`: The **`.gfa`** file used by GraphAligner.

- **Writable scratch/output location**
  - `scratch_path`: A directory path that **must already exist** and be writable.

- **Run parameters**
  - `flowcell`, `kit`: Guppy basecalling settings specifying the nanopore **flowcell type** and **sequencing kit** used for the run.
  - `fast5_chunksize`: Number of FAST5 files to group into each chunk for parallel basecalling.
  - `read_setsize`: Number of reads per “readset” FASTQ produced after filtering; used to split the full FASTQ stream into manageable batches for downstream alignment/barcode calling.
  - `min_readlength`, `max_readlength`: Read-length bounds used to filter basecalled reads before barcode calling.
  - `chunk_size`: Number of barcode groups to bundle for parallel per-barcode consensus calling (Medaka) and alignment.
  - `depth_threshold`: Minimum number of reads required for a barcode to be retained for grouping/consensus (filters out low-support barcodes).
  - `n_subsamples`: Number of per-barcode subsampling replicates to generate; should usually be set to `0` to indicate “no subsampling / use all reads”.

---

## What the Snakemake pipeline does

The `lDE26_Sequencing.smk` workflow takes raw nanopore **FAST5** and produces a final merged table (`output.tsv`) that, for each observed barcode group, includes:
- The called **barcode** (as a bitstring),
- The **per-barcode consensus sequence** (Medaka),
- The corresponding **per-barcode reference sequence** (derived from the GFA traversal), and alignment metadata (start position + CIGAR derived from minimap2 SAMs).

Conceptually, the pipeline proceeds in these stages:

1. **Chunk FAST5**  
   `scripts/group_fast5.py` splits the input FAST5 directory into `fast5_chunks/` for parallel basecalling.

2. **Basecall with Guppy (GPU)**  
   Each chunk is basecalled with `guppy_basecaller` into `fastq_chunks/` (configured via `flowcell` + `kit`).

3. **Create filtered readsets** *(checkpoint)*  
   `scripts/regroup_reads.py` filters reads by length (`min_readlength` / `max_readlength`) and writes fixed-size `reads/readset_*.fastq` files.

4. **Align readsets to the GFA graph**  
   GraphAligner maps reads to `ref.gfa`, producing `graph_output/readset_*.gaf`.

5. **Extract barcodes from alignments**  
   `scripts/get_barcodes.py` parses GAF output and writes per-read barcode calls (`graph_output/readset_*.tsv`).

6. **Choose barcode depth threshold + build codebook**  
   `scripts/set_barcode_threshold.py` counts barcode occurrences, saves a histogram plot, and writes an “inverse codebook” mapping `barcodeid → read list (+ barcode sequence)`.

7. **Group reads by barcode and write per-barcode references** *(checkpoint)*  
   `scripts/group_barcodes.py`:
   - Writes per-barcode FASTQs (`readgroups/`)
   - Writes per-barcode reference FASTAs (`grouprefs/`)
   - Optionally generates multiple **subsamples** (`subsample=0` means “use all reads”; positive values mean “randomly subsample N reads per barcode”)
   - Splits barcode groups into `chunk_*` folders for parallel consensus calling

8. **Call per-barcode consensus + align consensus back to reference**  
   For each per-barcode FASTQ:
   - `medaka_consensus` produces a consensus FASTA
   - `minimap2` aligns the consensus to the per-barcode reference, producing a SAM for downstream CIGAR parsing

9. **Aggregate per-subsample outputs**  
   The workflow concatenates consensus/reference FASTAs across chunks, aggregates SAM CIGARs, and writes:
   - `output_subsample={N}.tsv` (one per subsample)
   - `output.tsv` (merged across all subsamples)

### QC + export notebook

`Barcode_QC_and_Export.ipynb` is an example downstream analysis for the `lDE26` run. It demonstrates how to:

- Import `output.tsv` and basic run metadata from the config
- Compute barcode QC summaries such as:
  - **bit frequencies** across observed barcodes
  - **nearest-neighbor Hamming distance** distributions
  - **barcode depth** (reads-per-barcode) distributions
- Generate/export a cleaned table suitable for downstream analysis (e.g., joining barcode calls to a library design table to build a barcode→variant mapping)

---
