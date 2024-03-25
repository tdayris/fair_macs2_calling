rule deeptools_bamcoverage:
    input:
        unpack(get_deeptools_bamcoverage_input),
    output:
        protected(
            "results/{species}.{build}.{release}.dna/Coverage/{sample}.bw",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_mac2_calling/deeptools/bamcoverage/{species}.{build}.{release}.dna/{sample}.log",
    benchmark:
        "benchmark/fair_mac2_calling/deeptools/bamcoverage/{species}.{build}.{release}.dna/{sample}.tsv"
    params:
        genome="{build}",
        effective_genome_size=lambda wildcards: get_effective_genome_size(
            wildcards, genomes
        ),
        read_length=lambda wildcards: get_read_length(wildcards, samples),
        extra=dlookup(
            dpath="params/fair_mac2_calling/deeptools/bamcoverage",
            within=config,
            default="--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --ignoreForNormalization X Y MT",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/bamcoverage"


rule deeptools_plotcoverage:
    input:
        unpack(get_deeptools_plotcoverage_input),
    output:
        plot=report(
            "results/{species}.{build}.{release}.dna/Graphs/PlotCoverage.png",
            caption="../report/deeptools_plotcoverage.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_coverage",
                "species": "{species}.{build}.{release}",
            },
        ),
        raw_counts=temp(
            "tmp/fair_mac2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.raw"
        ),
        metrics=temp(
            "tmp/fair_mac2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.metrics"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * (60 * 5),
        tmpdir=tmp,
    log:
        "logs/fair_mac2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_mac2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna.tsv"
    params:
        dlookup(
            dpath="params/fair_mac2_calling/deeptools/plot_coverage",
            within=config,
            default="--skipZeros --coverageThresholds 1 --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/plotcoverage"


rule deeptools_fingerprint:
    input:
        unpack(get_deeptools_fingerprint_input),
    output:
        fingerprint=report(
            "results/{species}.{build}.{release}.dna/Graphs/PlotFingerprint.png",
            caption="../report/deeptools_plotfingerprint.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_fingerprint",
                "species": "{species}.{build}.{release}",
            },
        ),
        counts=temp(
            "tmp/fair_mac2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna/raw_counts.tab"
        ),
        qc_metrics=temp(
            "tmp/fair_mac2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna/qc_metrics.txt"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_mac2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_mac2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna.tsv"
    params:
        dlookup(
            dpath="params/fair_mac2_calling/deeptools/plot_fingerprint",
            within=config,
            default="--skipZeros --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/plotfingerprint"
