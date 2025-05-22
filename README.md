# !! Warning !!

This version of the pipeline has not been as thoroughly tested. There will be a new version coming soon that will be much more reproducible. This code is provided as-is to reflect the code used to generate the original results.

# SPIDR Paper Pipeline

This repository contains the computational pipeline used in the SPIDR (Split and Pool Identification of RBP targets) paper. SPIDR is a massively multiplexed method designed to simultaneously profile global RNA binding sites of dozens of RBPs in a single experiment.

## Overview

The pipeline processes SPIDR sequencing data through several key steps:
- Barcode identification and processing
- RNA alignment and analysis
- Cluster generation and analysis
- Quality control and visualization

## Requirements

- Snakemake
- Conda/Mamba
- Java (for barcode identification)
- Python dependencies (specified in environment files)

## Usage

The pipeline is implemented using Snakemake for workflow management. To run:

```bash
./run_pipeline.sh
```

This will execute the pipeline using the configuration specified in `config.yaml` and cluster settings in `cluster.yaml`.

## Pipeline Structure

- `scripts/`: Contains Python and Java processing scripts
- `envs/`: Conda environment specifications
- `Snakefile`: Main workflow definition
- `config.yaml`: Pipeline configuration
- `cluster.yaml`: Cluster configuration for distributed processing

## Citation

Please cite the following paper when using SPIDR:

```bibtex
@article {Wolin2023,
	author = {Erica Wolin and Jimmy K. Guo and Mario R. Blanco and Andrew A. Perez and Isabel N. Goronzy and Ahmed A. Abdou and Darvesh Gorhe and Mitchell Guttman and Marko Jovanovic},
	title = {SPIDR: a highly multiplexed method for mapping RNA-protein interactions uncovers a potential mechanism for selective translational suppression upon cellular stress},
	elocation-id = {2023.06.05.543769},
	year = {2023},
	doi = {10.1101/2023.06.05.543769},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {RNA binding proteins (RBPs) play crucial roles in regulating every stage of the mRNA life cycle and mediating non-coding RNA functions. Despite their importance, the specific roles of most RBPs remain unexplored because we do not know what specific RNAs most RBPs bind. Current methods, such as crosslinking and immunoprecipitation followed by sequencing (CLIP-seq), have expanded our knowledge of RBP-RNA interactions but are generally limited by their ability to map only one RBP at a time. To address this limitation, we developed SPIDR (Split and Pool Identification of RBP targets), a massively multiplexed method to simultaneously profile global RNA binding sites of dozens to hundreds of RBPs in a single experiment. SPIDR employs split-pool barcoding coupled with antibody-bead barcoding to increase the throughput of current CLIP methods by two orders of magnitude. SPIDR reliably identifies precise, single-nucleotide RNA binding sites for diverse classes of RBPs simultaneously. Using SPIDR, we explored changes in RBP binding upon mTOR inhibition and identified that 4EBP1 acts as a dynamic RBP that selectively binds to 5'-untranslated regions of specific translationally repressed mRNAs only upon mTOR inhibition. This observation provides a potential mechanism to explain the specificity of translational regulation controlled by mTOR signaling. SPIDR has the potential to revolutionize our understanding of RNA biology and both transcriptional and post-transcriptional gene regulation by enabling rapid, de novo discovery of RNA-protein interactions at an unprecedented scale.},
	URL = {https://www.biorxiv.org/content/early/2023/06/07/2023.06.05.543769},
	eprint = {https://www.biorxiv.org/content/early/2023/06/07/2023.06.05.543769.full.pdf},
	journal = {bioRxiv}
}
```

*Note: This citation will be updated when the final publication becomes available.* 