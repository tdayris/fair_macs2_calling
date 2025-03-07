rule fair_macs2_calling_datavzrd_homer_yaml:
    input:
        summary="tmp/fair_macs2_calling_summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv",
    output:
        yaml=temp(
            "tmp/fair_macs2_calling_datavzrd_homer_yaml/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.yaml"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_datavzrd_homer_yaml/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.config.log",
    benchmark:
        "benchmark/fair_macs2_calling_datavzrd_homer_yaml/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.config.tsv"
    params:
        sample="{sample}",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_datavzrd_yaml.py"


rule fair_macs2_calling_datavzrd_homer_render:
    input:
        config="tmp/fair_macs2_calling_datavzrd_homer_yaml/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.yaml",
        table="tmp/fair_macs2_calling_summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv",
    output:
        report(
            directory(
                "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}_reports"
            ),
            htmlindex="index.html",
            caption="../report/datavzrd_homer.rst",
            category="Peak Calling",
            subcategory="{macs2_peak_type}",
            labels={
                "species": "{species}.{build}.{release}",
                "sample": "{sample}",
            },
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_datavzrd_homer_render/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.render.log",
    benchmark:
        "benchmark/fair_macs2_calling_datavzrd_homer_render/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.render.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_datavzrd_homer_render", default=""
        ),
    wrapper:
        "v5.8.3/utils/datavzrd"
