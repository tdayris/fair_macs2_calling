rule xsv_cat:
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
        "v3.2.0/utils/xsv"


rule xsv_sort:
    input:
        table="tmp/csv/cat_rows/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp(
            "tmp/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv"
        ),
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
        "v3.2.0/utils/xsv"


rule xsv_fmt:
    input:
        table="tmp/csv/sort/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.csv",
    output:
        temp("tmp/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed"),
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
        "v3.2.0/utils/xsv"


rule bedtools_merge:
    input:
        "tmp/csv/fmt/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.bed",
    output:
        temp(
            "tmp/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed"
        ),
    log:
        "logs/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {}).get("bedtools", {}).get("merge_peaks", ""),
    wrapper:
        "v3.2.0/bio/bedtools/merge"
