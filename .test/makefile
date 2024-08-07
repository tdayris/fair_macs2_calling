SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt report.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/bedtools_merge_peaks.smk" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/datavzrd.smk" \
	"../workflow/rules/deeptools_bam_qc.smk" \
	"../workflow/rules/deeptools_correlations.smk" \
	"../workflow/rules/fair_fastqc_multiqc_pipeline.smk" \
	"../workflow/rules/fair_genome_indexer_pipeline.smk" \
	"../workflow/rules/fair_bowtie2_mapping_pipeline.smk" \
	"../workflow/rules/homer.smk" \
	"../workflow/rules/in_house_homer_post_process.smk" \
	"../workflow/rules/macs2.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/Snakefile"


# List of python files to lint
python_scripts := \
	"../workflow/scripts/build_datavzrd_yaml.py" \
	"../workflow/scripts/deeptools_multi_bigwig_summary_wrapper.py" \
	"../workflow/scripts/deeptools_plot_correlation_wrapper.py" \
	"../workflow/scripts/deeptools_plot_enrichment_wrapper.py" \
	"../workflow/scripts/deeptools_plot_pca_wrapper.py" \
	"../workflow/scripts/homer_to_json.py" \
	"../workflow/scripts/homer_wrapper.py" \
	"../workflow/scripts/merge_homer_summaries.py" \
	"../workflow/scripts/plot_homer_multisample.py"



.PHONY: all
all: report.txt


dag.png: format.txt linter.txt
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	snakemake -s ../workflow/Snakefile \
		--report report.zip 2>&1 | tee report.txt


pipeline.txt: linter.txt
	snakemake -s ../workflow/Snakefile \
		--cores 7 \
		--software-deployment-method conda \
		2>&1 | tee pipeline.txt


format.txt:
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt: black.txt format.txt
	snakemake -s ../workflow/Snakefile \
		--lint 2>&1 | tee linter.txt


black.txt:
	black $(python_scripts) 2>&1 | tee black.txt


.PHONY: clean
clean:
	rm -rvf linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt \
		   black.txt \
		   report.zip \
		   report
