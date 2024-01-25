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

| Step                                                       | Commands/Wrapper                                                                                                             |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| Download DNA Fasta from Ensembl                            | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-sequence.html)     |
| Remove non-canonical chromosomes                           | [pyfaidx](https://github.com/mdshw5/pyfaidx)                                                                         |
| Index DNA sequence                                         | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                         |
| Creatse sequence Dictionary                                | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html)          |
| Download GTF annotation                                    | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-annotation.html) |
| Fix format errors                                          | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                     |
| Remove non-canonical chromosomes, based on above DNA Fasta | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sq_filter_feature_from_fasta.html)                           |
| Remove `<NA>` Transcript support levels                    | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sp_filter_feature_by_attribute_value.html)                   |
| Download blacklisted regions                               | [Github source](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)                                            |
| Merge overlapping intervals                                | [bedtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bedtools/merge.html)                         |


### Raw-sequences QC with [`fair_fastqc_multiqc`](https://github.com/tdayris/fair_fastqc_multiqc/)

| Step    | Wrapper                                                                                      |
| ------- | -------------------------------------------------------------------------------------------- |
| FastQC  | [fastqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/fastqc.html)   |
| MultiQC | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html) |


### Bowtie2 Mapping with [`fair_bowtie2_mapping`](https://github.com/tdayris/fair_bowtie2_mapping/tree/main)

| Step             | Meta-Wrapper                                                                                                             | Wrapper                                                                                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| Bowtie2-build    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-build](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bowtie2/build.html)                                 |
| Fastp            |                                                                                                                          | [fastp](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastp.html)                                                 |
| Bowtie2-align    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [bowtie2-align](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/bowtie2/align.html)                                 |
| Sambamba sort    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-sort](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/sort.html)                                 |
| Sambamba-view    | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-view](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/view.html)                                 |
| Sambamba-markdup | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-markdup](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/markdup.html)                           |
| Sambamba-index   | [bowtie2-sambamba meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/bowtie2_sambamba.html) | [sambamba-index](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/sambamba/index.html)                               |
| Picard           |                                                                                                                          | [picard-collectmultiplemetrics](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html) |
| Samtools         |                                                                                                                          | [samtools-stats](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/stats.html)                               |
|


### Peak Calling

| Step              | Wrapper                                                                                                          |
| ----------------- | ---------------------------------------------------------------------------------------------------------------- |
| Coverage analysis | [deeptools-bamcoverage](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/bamcoverage.html) |
| Peak-Calling      | [macs2-callpeak](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/macs2/callpeak.html)               |
| Annotation        | [homer-annotatepeaks](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/homer/annotatePeaks.html)     |


### Quality Controls

| Step            | Wrapper/Official documentation                                                                                           |
| --------------- | ------------------------------------------------------------------------------------------------------------------------ |
| MultiQC         | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html)                             |
| Finger prints   | [deeptools-plotfingerprint](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/deeptools/plotfingerprint.html) |
| Plot Coverage   | [deeptools-plotcoverage](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/deeptools/plotcoverage.html)       |
| Plot Enrichment | [deeptools-plotenrichment](https://deeptools.readthedocs.io/en/develop/content/tools/plotEnrichment.html)                |
| Plot PCA        | [deeptools-plotpca](https://deeptools.readthedocs.io/en/develop/content/tools/plotPCA.html)                              |
