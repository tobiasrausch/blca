# Nanopore sequencing of bladder cancer

This repository contains the analysis code for the bioRxiv preprint [10.1101/2025.07.30.667694v1](https://www.biorxiv.org/content/10.1101/2025.07.30.667694v1).

Associated sequencing data have been deposited at the European Genome-phenome Archive (EGA) under the accession number EGAS00001008321.

## Installation

To create a mamba environment with the required tools just use

`git clone https://github.com/tobiasrausch/blca`

`cd blca/`

`make all`

## Downloading and preparing reference files

The genome files are not included in the repository but can be downloaded using

`cd genome/ && ./prepare_genome.sh`

## Basecalling and alignment

The alignments for the ONT data were generated using [Dorado](https://github.com/nanoporetech/dorado) version 0.6.0 using the most accurate PromethION basecalling model (sup) with modified base detection of 5-methylcytosine (5mC) and 5-hydroxymethylcytosine (5hmC) in CG contexts (sup,5mCG_5hmCG model). The base called data was then aligned to the human reference genome (GRCh38) using minimap2. These alignments were deposited at EGA (EGAS00001008321) and are assumed to be downloaded in the `alignment` subdirectory of this repository.

## Alignment statistics

To calculate the alignment error rate, genome coverage and other QC statistics.

`cd qc/ && ./qc.sh`

