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

As stated above, tumor-only sequencing approaches overestimate the TMB but short-reads have severe limitations in detecting SVs compared to long-reads, in particular for insertions. We therefore use an approach that implants the insertions into the reference to identify remaining rare and ultra-rare germline SVs in the matched short-read control data.

`cd consensusSV/genotypingSR/ && make all && ./genotypingSR.sh`

Afterwards you can recompute the consensus, including the additional filtering against the matched controls.

## Discovery of L1 fragments inserted at structural variant breakpoints

To discovery L1 fragments inserted at SV breakpoints we used [BreakTracer](https://github.com/tobiasrausch/breaktracer/releases).

`cd breaktracer/ && ./bt.sh`

## Mobile element insertion annotation

We used [SVAN](https://github.com/REPBIO-LAB/SVAN) and [Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF) to annotate all mobile element insertions.

`cd mei/ && make all && ./mei.sh && ./postprocess_SVAN.sh`

## ecDNA discovery

We employed [CoRAL](https://github.com/AmpliconSuite/CoRAL) to reconstruct ecDNA structures from long-read tumor data.

`cd ecDNA/ && make all && ./ecDNA.sh`

## Methylation calling

[Modkit](https://github.com/nanoporetech/modkit) was used to create pileup files of 5mC modified bases.

`cd modkit/ && ./modkit.sh`

## Short-read control data analysis

For the short-read control data, we used [Delly](https://github.com/dellytools/delly) to call germline structural variants, [MELT](https://melt.igs.umaryland.edu/) to discover germline mobile element insertions and the nf-core pipeline [sarek](https://github.com/nf-core/sarek) and [oncoanalyser](https://github.com/nf-core/oncoanalyser) for short variants (SNVs and InDels).

`cd sr-analysis/ && ./delly.sh`

`cd sr-analysis/ && ./melt.sh`

`cd sr-analysis/ && ./sarek.sh`

## PCAWG analysis overlapping CNV breakpoints with and without SVs with somatic L1 insertions

Using the short-read [PCAWG data](https://www.nature.com/articles/s41586-020-1969-6), we classified each copy-number variant (CNV) breakpoint as either SV-explained or SV-unexplained using a small overlap window of 5kbp. We also assessed the overlap between somatic L1 insertions for both categories, CNVs with and without SVs. The fraction of somatic L1s is significantly greater for CNV breakpoints without an SV, indicating a hidden landscape of SVs that involves L1 fragments at the SV breakpoint and these complex somatic SVs are apparently largely inaccessible to short-reads. For long-reads, [BreakTracer](https://github.com/tobiasrausch/breaktracer) resolves these complex SVs involving L1 sequence fragments.

`cd cnv_sv_l1/ && ./cnv.sh`

## Short-read bulk RNA-Seq

The short-read RNA-Seq data was analyzed with the [nf-core/rnaseq](https://github.com/nf-core/rnaseq/) pipeline to quantify the RNA gene expression and discover gene fusion candidates.

`cd rna-analysis/ && ./rnaseq.sh`

To quantify locus-specific L1 expression, we used the [L1EM](https://github.com/FenyoLab/L1EM) tool.

`cd rna-analysis/ && ./l1em.sh`

License
-------
All scripts are distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/blca/blob/main/LICENSE) file for more details.

Credits
-------
[HTSlib](https://github.com/samtools/htslib), [Samtools](https://github.com/samtools/samtools) and [Bcftools](https://github.com/samtools/bcftools) are heavily used for all genomic alignment and variant processing. [nf-core](https://github.com/nf-core) pipelines are employed for short-read variant calling and RNA-Seq data analysis. Most importantly, we thank all patients for participating in the study. 
