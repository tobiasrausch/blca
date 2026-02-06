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

## Read-depth profiles and copy-number variants

Genome-wide read-depth profiles

`cd coverage/ && ./coverage.sh`

## Single-nucleotide variant (SNVs) calling and small insertions and deletions (InDels)

SNVs and InDel calling using [ClairS](https://github.com/HKU-BAL/ClairS) and [DeepSomatic](https://github.com/google/deepsomatic).

`cd deepsomatic/ && ./deepsomatic.sh`

`cd clairS/ && make all && ./clairS.sh`

To compute the consensus of ClairS and DeepSomatic:

`cd consensusSNV && ./consensus.sh`

Tumor-only sequencing approaches overestimate the TMB (tumor mutational burden) due to rare and ultra-rare germline variants that are missing in panel-of-normals. Therefore, you can optionally genotype all variants using matched short-read control data.

`cd consensusSNV/genotypeSR/ && make all && ./genotypeSR.sh`

Afterwards you can recompute the consensus, including the additional filtering against the matched controls.

`cd consensusSNV && ./consensus.sh`

## Structural variant calling

For the ONT long-read data we applied [Severus](https://github.com/KolmogorovLab/Severus) and [Delly](https://github.com/dellytools/delly) to call structural variants with the [Schloissnig et al. data](https://www.nature.com/articles/s41586-025-09290-7) as a panel-of-normals.

`cd sv/ && ./sv.sh`

To compute the consensus of the Severus and Delly SV calls:

`cd consensusSV/ && ./somaticSV.sh`

## Discovery of L1 fragments inserted at structural variant breakpoints

To discovery L1 fragments inserted at SV breakpoints we used [BreakTracer](https://github.com/tobiasrausch/breaktracer/releases).

`cd breaktracer/ && ./bt.sh`

## Mobile element insertion annotation

We used [SVAN](https://github.com/REPBIO-LAB/SVAN) and [Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF) to annotate all mobile element insertions.

`cd mei/ && ./mei.sh`

## Methylation calling

[Modkit](https://github.com/nanoporetech/modkit) was used to create pileup files of 5mC modified bases.

`cd modkit/ && ./modkit.sh`
