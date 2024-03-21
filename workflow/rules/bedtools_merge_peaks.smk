rule fair_mac2_calling_xsv_cat_macs2_peaks:
    input:
        table=expand(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed",
            sample=samples.sample_id,
            allow_missing=True,
        ),
    output:
        temp(
            "tmp/fair_mac2_calling/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_mac2_calling/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_mac2_calling/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="cat rows",
        extra=dlookup(dpath="params/fair_macs2_calling/xsv/cat_macs2_peaks", within=config, default="--no-headers --delimiter $'\t'"),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


rule fair_mac2_calling_xsv_sort_macs2_concat_peaks:
    input:
        table="tmp/fair_mac2_calling/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp(
            "tmp/fair_mac2_calling/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_mac2_calling/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_mac2_calling/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="sort",
        extra=dlookup(dpath="params/fair_macs2_calling/xsv/sort_macs2_peaks", within=config, default="--no-headers --numeric --select 1-3"),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


rule fair_mac2_calling_xsv_fmt_macs2_sorted_peaks:
    input:
        table="tmp/fair_mac2_calling/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp("tmp/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="fmt",
        extra=dlookup(dpath="params/fair_macs2_calling/xsv/format_macs2_peaks", within=config, default="--out-delimiter $'\t'"),
    wrapper:
        f"{snakemake_wrappers_prefix}/utils/xsv"


rule fair_mac2_calling_bedtools_sort_macs2_concat_beds:
    input:
        in_file="tmp/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed",
    output:
        temp(
            "tmp/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.sorted.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 8),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/fair_mac2_calling/bdtools/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_mac2_calling/bedtools/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(dpath="params/fair_macs2_calling/bedtools/sort_macs2_bed", within=config, default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bedtools/sort"


rule fair_mac2_calling_bedtools_merge_macs2_sorted_peaks:
    input:
        "tmp/fair_mac2_calling/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.sorted.bed",
    output:
        temp(
            "tmp/fair_mac2_calling/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/fair_mac2_calling/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_mac2_calling/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(dpath="params/fair_macs2_calling/bedtools/merge_macs2_peaks", within=config, default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/bedtools/merge"
