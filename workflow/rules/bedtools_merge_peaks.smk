rule xsv_cat_macs2_peaks:
    input:
        table=expand(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed",
            sample=samples.sample_id,
            allow_missing=True,
        ),
    output:
        temp(
            "tmp/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="cat rows",
        extra="--no-headers --delimiter $'\t'",
    wrapper:
        "v3.3.3/utils/xsv"


rule xsv_sort_macs2_concat_peaks:
    input:
        table="tmp/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp(
            "tmp/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="sort",
        extra="--no-headers --numeric --select 1-3",
    wrapper:
        "v3.3.3/utils/xsv"


rule xsv_fmt_macs2_sorted_peaks:
    input:
        table="tmp/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp("tmp/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
    group:
        "concat_beds"
    log:
        "logs/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        subcommand="fmt",
        extra="--out-delimiter $'\t'",
    wrapper:
        "v3.3.3/utils/xsv"


rule bedtools_merge_macs2_sorted_peaks:
    input:
        "tmp/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed",
    output:
        temp(
            "tmp/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 4),
        runtime=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {}).get("bedtools", {}).get("merge_peaks", ""),
    wrapper:
        "v3.3.3/bio/bedtools/merge"
