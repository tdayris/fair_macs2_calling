# 1.0.2

## Features

* Missing parameters added to documentation

## Fixes

* Typo in usage
* Fix csv.Sniffer too large data content

# 1.0.1

## Features

* Include FastQC in quality report
* Snakemake v8+ compatibility
* Snakemake-wrappers v3.3.3
* All MultiQC are in same tab in the report
* Snakemake environment always point to the latest version of Snakemake

## Fixes

* Typo in multiqc tab of snakemake report
* Fix x/y axis visibility in catplot



# 1.0.0

## Features

* deeptools used to compute coverage and perform QC graphs
* macs2 calls peaks in both narrow and broad modes
* homer annotates peaks
* deeptools plotcoverage and plotfingerprint
* integration of fair_bowtie2_mapping
* integration of fair_genome_indexer
* Snakemake workflow complete
* deeptools correlation and PCA over all peaks
* in-house scripts + datavzrd for Homer results