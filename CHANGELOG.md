# 2.0.2

## Features

* Prevent homer temp file collision
* Update github checkout images

## Fix

* Missing read length in default genome size computation

# 2.0.1

## Features

* Resources reservation

## Fixes

* Samtools configuration

# 2.0.0

Breaking change: Non-canonical chromosomes removed by defaults

## Features:

* fair_genome_indexer update to [3.0.0](https://github.com/tdayris/fair_genome_indexer/releases/tag/3.0.0)
* fair_bowtie2_mapping update to [3.0.0](https://github.com/tdayris/fair_bowtie2_mapping/releases/tag/3.0.0)
* snakemake-wrappers updated to [v3.3.3](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/index.html)
* snakemake github action updated to [1.25.1](https://github.com/snakemake/snakemake-github-action/releases/tag/v1.25.1)

## Fix:

* Wrong parameters passed to MultiQC when importing pipeline with no parameters
* Report testing deactivated as long as TBD issue is opened at snakeamke

## Documentation:

* Pipeline description updated
* Usage generalized
* Gustave Roussy users have a dedicated usage section

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