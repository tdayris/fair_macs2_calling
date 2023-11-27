rule deeptools_bamcoverage:
    input:
        bam="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam",
        bai="results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai",
        blacklist="reference/blacklist/{species}.{build}.{release}.merged.bed",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw"
        )
    log:
        "logs/deeptools/bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.log"
    benchmark:
        "benchmark/deeptools/bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        genome="{build}",
        read_length="100",
        extra=config.get("params", {}).get("deeptools", {}).get("bamcoverage", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/bamcoverage"


rule deeptools_plotcoverage:
    input:
        unpack(get_deeptools_plotcoverage_input)
    output:
        plot="results/{species}.{build}.{release}.{datatype}/PlotCoverage/{sample}.png",
        raw_counts=temp("tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}_coverage.raw"),
        metrics=temp("tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}_coverage.metrics"),
    log:
        "logs/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}.log"
    benchmark:
        "benchmark/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        extra=config.get("params", {}).get("deeptools", {}).get("plotcoverage", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/deeptools/plotcoverage"