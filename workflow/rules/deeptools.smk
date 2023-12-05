rule deeptools_bamcoverage:
    input:
        unpack(get_deeptools_bamcoverage_input),
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
        ),
    log:
        "logs/deeptools/bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/deeptools/bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        genome="{build}",
        effective_genome_size=lambda wildcards: get_effective_genome_size(
            wildcards, genomes
        ),
        extra=config.get("params", {}).get("deeptools", {}).get("bamcoverage", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/bamcoverage"


rule deeptools_plotcoverage:
    input:
        unpack(get_deeptools_plotcoverage_input),
    output:
        plot=report(
            "results/{species}.{build}.{release}.{datatype}/PlotCoverage.png",
            caption="../report/deeptools_plotcoverage.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_coverage",
                "species": "{species}.{build}.{release}",
            },
        ),
        raw_counts=temp(
            "tmp/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}/Coverage.raw"
        ),
        metrics=temp(
            "tmp/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}/Coverage.metrics"
        ),
    log:
        "logs/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=config.get("params", {}).get("deeptools", {}).get("plot_coverage", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotcoverage"


rule deeptools_fingerprint:
    input:
        unpack(get_deeptools_fingerprint_input),
    output:
        plot=report(
            "results/{species}.{build}.{release}.{datatype}/PlotFingerprint.png",
            caption="../report/deeptools_plotfingerprint.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_fingerprint",
                "species": "{species}.{build}.{release}",
            },
        ),
        raw_counts=temp(
            "tmp/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}/raw_counts.tab"
        ),
        metrics=temp(
            "tmp/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}/qc_metrics.txt"
        ),
    log:
        "logs/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}.tsv"
    params:
        config.get("params", {}).get("deeptools", {}).get("plot_fingerprint", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotfingerprint"
