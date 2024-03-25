rule macs2_callpeak_narrow:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_peaks.xls"
        ),
        temp(
            ensure(
                "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_peaks.narrowPeak",
                non_empty=True,
            )
        ),
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_summits.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{sample}.narrow.log",
    benchmark:
        "benchmark/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{sample}.narrow.tsv"
    params:
        lambda wildcards: get_macs2_callpeak_params(wildcards, samples, config),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/macs2/callpeak"


use rule macs2_callpeak_narrow as macs2_callpeak_broad with:
    output:
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.xls"
        ),
        temp(
            ensure(
                "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.broadPeak",
                non_empty=True,
            )
        ),
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.gappedPeak"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{sample}.broad.log",
    benchmark:
        "benchmark/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{sample}.broad.tsv"


rule macs2_peaks_to_csv:
    input:
        table="tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}_peaks.{macs2_peak_type}",
    output:
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.csv"
        ),
    log:
        "logs/fair_macs2_calling/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    group:
        "macs2_reformat"
    params:
        subcommand="select",
        extra=dlookup(dpath="params/fair_macs2_calling/xsv/macs_peaks_to_csv", within=config, default="1-6 --delimiter $'\t'"),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


rule macs2_csv_to_bed:
    input:
        table="tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.csv",
    output:
        temp(
            "tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    group:
        "macs2_reformat"
    params:
        subcommand="fmt",
        extra=dlookup(dpath="params/fair_macs2_calling/xsv/macs_csv_to_bed", within=config, default="--out-delimiter $'\t'"),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


rule sort_macs2_bed:
    input:
        in_file="tmp/fair_macs2_calling/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.bed",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/bdtools/sort/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/bedtools/sort/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    params:
        extra=dlookup(dpath="params/fair_macs2_calling/bedtools/sort", within=config, default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bedtools/sort"
