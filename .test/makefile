SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt report.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/deeptools.smk" \
	"../workflow/rules/homer.smk" \
	"../workflow/rules/macs2.smk" \
	"../workflow/rules/fair_genome_indexer.smk" \
	"../workflow/rules/fair_bowtie2_mapping.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/Snakefile"


.PHONY: all
all: report.txt


dag.png: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--report report.html 2>&1 | tee report.txt


pipeline.txt: format.txt linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--cores 2 \
		--use-conda \
		--use-singularity \
		--workflow-profile ../profiles 2>&1 | tee pipeline.txt


format.txt:
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt:
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--lint 2>&1 | tee linter.txt


.PHONY: clean
clean:
	rm -rv linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt