This pipeline requires two configuration file:

# `config.yaml`

A standard `Snakemake` configuration, yaml-formatted file containing a list of
all parameters accepted in this workflow:

* `samples`: Path to the file containing link between samples and their fastq file(s)
* `params`: Per-tool list of optional parameters

Example:

```
samples: config/samples.csv

# Optional parameters
params:
  # Optional parameters for FastQC
  fastqc: ""
  # Optional parameters for FastP
  fastp:
    # Optional command line adapters
    adapters: ""
    # Optional command line parameters
    extra: ""
  bowtie2:
    # Optional parameters for bowtie2-build
    build: ""
    # Optional parameters for bowtie2-align
    align: ""
  sambamba:
    # Optional parameters for sambamba view
    view: "--format 'bam' --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' "
    # Optional parameters for sambamba markdup
    markdup: "--remove-duplicates --overflow-size 500000"
  picard:
    # Mapping QC optional parameters
    metrics: ""
  # Optional parameters for samtools stats
  samtools: ""
  macs2:
    # Optional parameters for Macs2 Callpeaks
    callpeak: ""
  homer:
    # Optional parameters for Homer AnnotatePeaks
    annotatepeaks: "-homer2"
  deeptools:
    # Optional parameters for DeepTools bamCoverage
    bamcoverage: "--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --ignoreForNormalization X Y MT"
    # Optional parameters for DeepTools plotCoverage
    plot_coverage: "--skipZeros --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4"
    # Optional parameters for DeepTools plotFingerprint
    plot_fingerprint: "--skipZeros --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4"
    # Optional parameters for DeepTools plotPCA
    plot_pca: "--ntop 1000"
    # Optional parameters for DeepTools plotEnrichment
    plot_enrichment: "--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --smartLabels"
    # Optional parameters for DeepTools plotCorrelations
    plot_correlation: "--whatToPlot heatmap --corMethod spearman --skipZeros --plotNumbers --colorMap RdYlBu"
  # Optional parameters for MultQC
  multiqc: "--module deeptools --module macs2 --module picard --module fastqc --module fastp --module samtools --module bowtie2 --module sambamba --zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"
```

# `samples.csv`

A CSV-formatted text file containing the following mandatory columns:

* sample_id: Unique name of the sample
* upstream_file: Path to upstream fastq file
* species: The species name, according to Ensembl standards
* build: The corresponding genome build, according to Ensembl standards
* release: The corresponding genome release, according to Ensembl standards
* downstream_file: Optional path to downstream fastq file
* input: Sample id of the corresponding input signal

Example:

```
sample_id,upstream_file,downstream_file,species,build,release,input
sac_a,data/reads/a.scerevisiae.1.fq,data/reads/a.scerevisiae.2.fq,saccharomyces_cerevisiae,R64-1-1,110,sac_a_input
sac_a_input,data/reads/a.scerevisiaeI.1.fq,data/reads/a.scerevisiaeI.2.fq,saccharomyces_cerevisiae,R64-1-1,110,
```

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.