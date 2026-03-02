## Overview
SpaceA-Paper is a repository containing the complete set of code, configuration files, and visualization results for the research paper centered on the SpaceA technology. SpaceA is a novel single-cell spatial chromatin conformation capture technology that enables high-resolution analysis of chromatin interactions in the spatial context of embryonic/tissue samples, revealing the spatial regulatory mechanisms of chromatin structure during development.

This repository organizes the analytical workflows and results by functional modules, covering data preprocessing, spot/cluster/read filtering, Higashi-based chromatin interaction analysis, and multi-dimensional visualization of results (e.g., chromatin contact distance, cell cycle correlation, UMAP spatial clustering).

## Repository Structure
### 1_Pre-processing
Core module for raw sequencing data preprocessing of spaSPRITE (the experimental basis of SpaceA) samples:
- `Snakefile`: Snakemake workflow configuration file that defines the end-to-end preprocessing pipeline (including raw read quality control, barcode parsing, alignment to reference genome, duplicate removal, etc.).
- `config.yaml`: Global parameter configuration for preprocessing (e.g., reference genome path, quality control thresholds, output directory).
- Sample-specific config files (`config_spaSPRITE_E11.5L1&E11.5L2.txt`, etc.): Parameter files tailored to different embryonic stage/tissue samples (E11.5L1, E12.5L5, E13.5C1, etc.), defining sample-specific alignment and filtering rules.
- `split_dpm_from_full_BC_fq.py`: Python script to split DPM (DNA Proximity Module) related reads from fastq files with full barcodes, a key step for parsing spaSPRITE sequencing data.

### 2_Filter_SpotsClustersReads
Module for quality filtering of spots, clusters, and sequencing reads:
- Implements filtering strategies for low-quality spots (e.g., spots with insufficient read depth), abnormal clusters (e.g., clusters with no biological significance), and invalid reads (e.g., unaligned/low-quality reads).
- Outputs filtered datasets for downstream Higashi analysis and visualization, ensuring the reliability of subsequent results.

### 3_Higashi
Module for chromatin interaction analysis using the FastHigashi algorithm:
- `FastHigashi/`: Directory storing FastHigashi source code dependencies, configuration files, and intermediate results (FastHigashi is an optimized Higashi tool for single-cell spatial chromatin interaction analysis).
- `higashi_pre.ipynb`: Jupyter Notebook for Higashi preprocessing, including data format conversion (adapting SpaceA data to Higashi input requirements), parameter tuning, and initial normalization of chromatin interaction matrices.

### Figs
Central directory for storing all visualization results of analytical findings, organized by analysis themes:
- `Cellcycle/`: Visualizations of cell cycle phase distribution and its correlation with chromatin interactions.
- `Contact_Distance/`: Analysis of chromatin contact distance distribution (global and tissue-specific).
- `Contact_Distance_Sample/`: Sample-specific chromatin contact distance analysis and cross-sample comparison.
- `Contact_Heatmap/`: Heatmaps of chromatin interaction matrices (SpaceA data-only).
- `Contact_Num/`: Statistics and visualization of chromatin interaction counts across samples/clusters.
- `Filter_Spot_Sprite/`: Visualization of spot filtering effects (before/after filtering comparison).
- `Heatmap_SpaceA_BulkHiC/`: Comparative heatmaps of chromatin interactions detected by SpaceA vs. Bulk Hi-C.
- `Lineplot_E1value/`: Line charts of E1 value (a quantitative index for chromatin interaction/expression) dynamics across developmental stages/clusters.
- `NoSwitchRange_Expression/`: Gene expression analysis in "NoSwitchRange" genomic regions and its correlation with chromatin interactions.
- `Pseudotime_Analysis/`: Pseudotime trajectory analysis of cells/spots, showing dynamic changes of chromatin interactions during development.
- `UMAP_EmbryoSpatial/`: UMAP dimensionality reduction visualization integrating spatial location and chromatin interaction features of embryonic samples.

### Root Directory Files
- `README-CH.md`: Chinese README file.
- `config.yaml`: Global configuration file for the entire repository, defining shared parameters across modules.

## Notes
- All scripts are optimized for SpaceA/spaSPRITE data of embryonic mouse samples; for other species/tissues, adjust reference genome paths, filtering thresholds, and algorithm parameters accordingly.
- Ensure consistent versioning of dependencies (e.g., FastHigashi, Snakemake) to avoid analytical errors caused by version incompatibility.
- Raw sequencing data used in this study is available from the corresponding NCBI/SRA database (accession number to be added in the paper).

## Citation
If you use the code or results from this repository in your research, please cite the SpaceA technology research paper (citation information to be updated upon publication).

## Contact
For technical issues or questions about the repository, please contact Yuetong XU (corresponding author) via email (email address to be provided in the paper).
