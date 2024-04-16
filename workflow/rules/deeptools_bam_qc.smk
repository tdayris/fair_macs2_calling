rule fair_macs2_calling_deeptools_bamcoverage:
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
        "logs/fair_macs2_calling/deeptools/bamcoverage/{species}.{build}.{release}.dna/{sample}.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/bamcoverage/{species}.{build}.{release}.dna/{sample}.tsv"
    params:
        genome="{build}",
        effective_genome_size=lambda wildcards: get_effective_genome_size(
            wildcards, genomes
        ),
        read_length=lambda wildcards: get_read_length(wildcards, samples),
        extra=lookup_config(
            dpath="params/fair_macs2_calling/deeptools/bamcoverage",
            default="--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --ignoreForNormalization X Y MT",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/bamcoverage"


rule fair_macs2_calling_deeptools_plotcoverage:
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
            "tmp/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.raw"
        ),
        metrics=temp(
            "tmp/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.metrics"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * (60 * 5),
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling/deeptools/plot_coverage",
            default="--skipZeros --coverageThresholds 1 --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/plotcoverage"


rule fair_macs2_calling_deeptools_fingerprint:
    input:
        bam_files=collect(
            "results/{sample.species}.{sample.build}.{sample.release}.dna/Mapping/{sample.sample_id}.bam",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bam_idx=collect(
            "results/{sample.species}.{sample.build}.{sample.release}.dna/Mapping/{sample.sample_id}.bam.bai",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
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
            "tmp/fair_macs2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna/raw_counts.tab"
        ),
        qc_metrics=temp(
            "tmp/fair_macs2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna/qc_metrics.txt"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna.tsv"
    params:
        lookup_config(
            dpath="params/fair_macs2_calling/deeptools/plot_fingerprint",
            default="--skipZeros --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/deeptools/plotfingerprint"
