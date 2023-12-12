[![Snakemake](https://img.shields.io/badge/snakemake-≥7.29.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_macs2_calling/workflows/Tests/badge.svg?branch=main)](https://github.com/tdayris/fair_macs2_calling/actions?query=branch%3Amain+workflow%3ATests)

Do not use. Active dev.

Snakemake workflow used to call peaks with Macs2 and annotate them with Homer.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_macs2_calling) 
it is also available [locally](https://github.com/tdayris/fair_macs2_calling/blob/main/workflow/report/usage.rst) on a single page.
 
## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_macs2_calling/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_macs2_calling/blob/main/workflow/report/material_methods.rst) textually. Web-links are available below:

![workflow_rulegraph](dag.png)

### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/tree/main)

| Step                          | Wrapper - Script                                                                                                                     |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| Download DNA fasta            | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/reference/ensembl-sequence.html)                     |
| Download cDNA fasta           | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/reference/ensembl-sequence.html)                     |
| Download GTF annotation       | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/reference/ensembl-annotation.html)                 |
| Samtools index fasta          | [samtools-faidx](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/samtools/faidx.html)                                   |
| Picard sequence dictionary    | [picard-createsequencedictionary](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/picard/createsequencedictionary.html) |
| Download VCF variation        | [ensembl-variation](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/reference/ensembl-variation.html)                   |
| Fix Ensembl GTF common errors | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                                     |
| Download known blacklist      | [Github source](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)                                                            |

### Bowtie2 Mapping with [`fair_bowtie2_mapping`](https://github.com/tdayris/fair_bowtie2_mapping/tree/main)

| Step             | Meta-Wrapper                                                                                                             | Wrapper                                                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| Bowtie2-build    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/bowtie2/build.html)                                 |
| Fastp            |                                                                                                                          | [fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastp.html)                                                 |
| Bowtie2-align    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/bowtie2/align.html)                                 |
| Sambamba sort    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/sambamba/sort.html)                                 |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/sambamba/view.html)                                 |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/sambamba/markdup.html)                           |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/sambamba/index.html)                               |
| Picard           |                                                                                                                          | [picard-collectmultiplemetrics](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html) |


### Peak Calling

| Step              | Wrapper                                                                                                          |
| ----------------- | ---------------------------------------------------------------------------------------------------------------- |
| Coverage analysis | [deeptools-bamcoverage](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/bamcoverage.html) |
| Peak-Calling      | [macs2-callpeak](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/macs2/callpeak.html)               |
| Annotation        | [homer-annotatepeaks](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/homer/annotatePeaks.html)     |


### Quality Controls

| Step          | Wrapper                                                                                                                  |
| ------------- | ------------------------------------------------------------------------------------------------------------------------ |
| MultiQC       | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/multiqc.html)                             |
| Finger prints | [deeptools-plotfingerprint](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/deeptools/plotfingerprint.html) |
| Plot Coverage | [deeptools-plotcoverage](https://snakemake-wrappers.readthedocs.io/en/v3.0.0/wrappers/deeptools/plotcoverage.html)       |
