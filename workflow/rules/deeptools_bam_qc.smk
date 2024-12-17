rule fair_macs2_calling_deeptools_bamcoverage:
    input:
        bam=lambda wildcards: get_aln(wildcards),
        bai=lambda wildcards: get_aln(wildcards, index=True),
        blacklist=lambda wildcards: get_blacklist(wildcards, genomes=genomes),
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling_deeptools_bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_bamcoverage/{species}.{build}.{release}.{datatype}/{sample}.tsv"
    params:
        genome="{build}",
        effective_genome_size=lambda wildcards: get_effective_genome_size(
            wildcards, genomes
        ),
        read_length=lambda wildcards: get_read_length(wildcards, samples),
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_bamcoverage",
            default="--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --ignoreForNormalization X Y MT",
        ),
    wrapper:
        "v5.5.0/bio/deeptools/bamcoverage"


rule fair_macs2_calling_deeptools_plotcoverage:
    input:
        unpack(get_deeptools_plotcoverage_input),
    output:
        plot=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/PlotCoverage.png",
            caption="../report/deeptools_plotcoverage.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_coverage",
                "species": "{species}.{build}.{release}.{datatype}",
            },
        ),
        raw_counts=temp(
            "tmp/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.raw"
        ),
        metrics=temp(
            "tmp/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.metrics"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * (60 * 5),
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_plotcoverage",
            default="--skipZeros --coverageThresholds 1 --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        "v5.5.0/bio/deeptools/plotcoverage"


rule fair_macs2_calling_deeptools_fingerprint:
    input:
        bam_files=collect(
            "results/{sample.species}.{sample.build}.{sample.release}.{datatype}/Mapping/{sample.sample_id}.bam",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        bam_idx=collect(
            "results/{sample.species}.{sample.build}.{sample.release}.{datatype}/Mapping/{sample.sample_id}.bam.bai",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
    output:
        fingerprint=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/PlotFingerprint.png",
            caption="../report/deeptools_plotfingerprint.rst",
            category="Coverage analysis",
            subcategory="Coverage",
            labels={
                "figure": "plot_fingerprint",
                "species": "{species}.{build}.{release}",
            },
        ),
        counts=temp(
            "tmp/fair_macs2_calling_deeptools_fingerprint/{species}.{build}.{release}.{datatype}/raw_counts.tab"
        ),
        qc_metrics=temp(
            "tmp/fair_macs2_calling_deeptools_fingerprint/{species}.{build}.{release}.{datatype}/qc_metrics.txt"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling_deeptools_fingerprint/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_fingerprint/{species}.{build}.{release}.{datatype}.tsv"
    params:
        lookup_config(
            dpath="params/fair_macs2_calling_deeptools_fingerprint",
            default="--skipZeros --ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4",
        ),
    wrapper:
        "v5.5.0/bio/deeptools/plotfingerprint"
