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
  bowtie2:
    # Optional parameters for bowtie2-build
    build: ""
    # Optional parameters for bowtie2-align
    align: ""
  sambamba:
    # Optional parameters for sambamba view
    view: "--format 'bam' --filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' "
    # Optional parameters for sambamba markdup
    markdup: "--remove-duplicates"
  picard:
    # Mapping QC optional parameters
    metrics: ""
  macs2:
    # Optional parameters for Macs2 Callpeaks
    callpeak: ""
  homer:
    # Optional parameters for Homer AnnotatePeaks
    annotatepeaks: ""
  deeptools:
    # Optional parameters for DeepTools bamCoverage
    bamcoverage: ""
    # Optional parameters for DeepTools plotCoverage
    plot_coverage: ""
    # Optional parameters for DeepTools plotFingerprint
    plot_fingerprint: ""
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