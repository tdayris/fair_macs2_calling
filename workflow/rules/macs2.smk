rule macs2_callpeak_narrow:
    input:
        unpack(get_macs2_callpeak_input),
    output:
        temp(
            "tmp/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_peaks.xls"
        ),
        temp(
            ensure(
                "tmp/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_peaks.narrowPeak",
                non_empty=True,
            )
        ),
        temp(
            "tmp/macs2/{species}.{build}.{release}.{datatype}/narrowPeak/{sample}_summits.bed"
        ),
    log:
        "logs/macs2/{species}.{build}.{release}.{datatype}/{sample}.narrow.log",
    benchmark:
        "benchmark/macs2/{species}.{build}.{release}.{datatype}/{sample}.narrow.tsv"
    params:
        lambda wildcards: get_macs2_callpeak_params(wildcards, samples, config),
    wrapper:
        "v3.3.3/bio/macs2/callpeak"


use rule macs2_callpeak_narrow as macs2_callpeak_broad with:
    output:
        temp(
            "tmp/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.xls"
        ),
        temp(
            ensure(
                "tmp/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.broadPeak",
                non_empty=True,
            )
        ),
        temp(
            "tmp/macs2/{species}.{build}.{release}.{datatype}/broadPeak/{sample}_peaks.gappedPeak"
        ),
    log:
        "logs/macs2/{species}.{build}.{release}.{datatype}/{sample}.broad.log",
    benchmark:
        "benchmark/macs2/{species}.{build}.{release}.{datatype}/{sample}.broad.tsv"


rule macs2_peaks_to_csv:
    input:
        table="tmp/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}_peaks.{macs2_peak_type}",
    output:
        temp(
            "tmp/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.csv"
        ),
    log:
        "logs/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    group:
        "macs2_reformat"
    params:
        subcommand="select",
        extra="1-6 --delimiter $'\t'",
    wrapper:
        "v3.3.3/utils/xsv"


rule macs2_csv_to_bed:
    input:
        table="tmp/macs2/{species}.{build}.{release}.{datatype}/{macs2_peak_type}/{sample}.csv",
    output:
        protected(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed"
        ),
    log:
        "logs/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/xsv/select/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    group:
        "macs2_reformat"
    params:
        subcommand="fmt",
        extra="--out-delimiter $'\t'",
    wrapper:
        "v3.3.3/utils/xsv"
