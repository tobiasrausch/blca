SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .coral
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y --override-channels -c conda-forge -c bioconda samtools bcftools bedtools htslib cnvkit pandas && pip install --force-reinstall statsmodels==0.14.0 pydeseq2==0.4.3 && touch .tools

.coral: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && pip install cvxopt intervaltree matplotlib numpy pandas pysam gurobipy==11.0.3 && touch .coral

clean:
	rm -rf mamba/ ${TARGETS}

