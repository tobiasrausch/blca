SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .clairSTO
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba create -n clairs-to -y --override-channels -c conda-forge -c pytorch -c bioconda pytorch tqdm clair3 bcftools einops scipy scikit-learn python=3.9.0 samtools bcftools bedtools htslib whatshap longphase && touch .tools

.clairSTO: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && source activate clairs-to && git clone https://github.com/HKU-BAL/ClairS-TO.git && mkdir -p mamba/envs/clairs-to/bin/clairs-to_models && mkdir -p mamba/envs/clairs-to/bin/clairs-to_databases && mkdir -p mamba/envs/clairs-to/bin/clairs-to_cna_data && wget http://www.bio8.cs.hku.hk/clairs-to/models/clairs-to_models.tar.gz && wget http://www.bio8.cs.hku.hk/clairs-to/databases/clairs-to_databases.tar.gz && wget http://www.bio8.cs.hku.hk/clairs-to/cna_data/reference_files.tar.gz && tar -zxvf clairs-to_models.tar.gz -C mamba/envs/clairs-to/bin/clairs-to_models/ && tar -zxvf clairs-to_databases.tar.gz -C mamba/envs/clairs-to/bin/clairs-to_databases/ && tar -zxvf reference_files.tar.gz -C mamba/envs/clairs-to/bin/clairs-to_cna_data/ && rm *.tar.gz && touch .clairSTO

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/
