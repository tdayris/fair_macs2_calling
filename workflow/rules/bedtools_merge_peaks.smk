rule fair_macs2_calling_xsv_cat_macs2_peaks:
    input:
        table=expand(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed",
            sample=samples.sample_id,
            allow_missing=True,
        ),
    output:
        temp(
            "tmp/fair_macs2_calling_xsv_cat_macs2_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_macs2_calling_xsv_cat_macs2_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_xsv_cat_macs2_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="cat rows",
        extra=lookup_config(
            dpath="params/fair_macs2_calling_xsv_cat_macs2_peaks",
            default="--no-headers --delimiter $'\t'",
        ),
    wrapper:
        "v5.8.3/utils/xsv"


rule fair_macs2_calling_xsv_sort_macs2_concat_peaks:
    input:
        table="tmp/fair_macs2_calling_xsv_cat_macs2_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp(
            "tmp/fair_macs2_calling_xsv_sort_macs2_concat_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_macs2_calling_xsv_sort_macs2_concat_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_xsv_sort_macs2_concat_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="sort",
        extra=lookup_config(
            dpath="params/fair_macs2_calling_xsv_sort_macs2_concat_peaks",
            default="--no-headers --numeric --select 1-3",
        ),
    wrapper:
        "v5.8.3/utils/xsv"


rule fair_macs2_calling_xsv_fmt_macs2_sorted_peaks:
    input:
        table="tmp/fair_macs2_calling_xsv_sort_macs2_concat_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp(
            "tmp/fair_macs2_calling_xsv_fmt_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/fair_macs2_calling_xsv_fmt_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_xsv_fmt_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="fmt",
        extra=lookup_config(
            dpath="params/fair_macs2_calling_xsv_fmt_macs2_sorted_peaks",
            default="--out-delimiter $'\t'",
        ),
    wrapper:
        "v5.8.3/utils/xsv"


rule fair_macs2_calling_bedtools_sort_macs2_concat_beds:
    input:
        in_file="tmp/fair_macs2_calling_xsv_fmt_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed",
    output:
        temp(
            "tmp/fair_macs2_calling_bedtools_sort_macs2_concat_beds/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.sorted.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 8),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_bedtools_sort_macs2_concat_beds/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_bedtools_sort_macs2_concat_beds/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_bedtools_sort_macs2_concat_beds",
            default="",
        ),
    wrapper:
        "v5.8.3/bio/bedtools/sort"


rule fair_macs2_calling_bedtools_merge_macs2_sorted_peaks:
    input:
        "tmp/fair_macs2_calling_bedtools_sort_macs2_concat_beds/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.sorted.bed",
    output:
        temp(
            "tmp/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks",
            default="",
        ),
    wrapper:
        "v5.8.3/bio/bedtools/merge"
