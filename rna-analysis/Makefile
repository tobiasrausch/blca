SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -n clairs-to -y --override-channels -c conda-forge -c bioconda python=2.7.15 bwa=0.7.17 samtools=1.9 numpy=1.14.3 scipy=1.1.0 pysam=0.15.0 bedtools=2.27.1 && touch .tools

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/
