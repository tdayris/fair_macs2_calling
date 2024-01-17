rule datavzrd_homer_yaml:
    input:
        summary="tmp/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv",
    output:
        yaml=temp(
            "tmp/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.yaml"
        ),
    log:
        "logs/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.config.log",
    benchmark:
        "benchmark/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.config.tsv"
    params:
        sample="{sample}",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_datavzrd_yaml.py"


rule datavzrd_homer_render:
    input:
        config="tmp/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.yaml",
        table="tmp/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv",
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
    log:
        "logs/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.render.log",
    benchmark:
        "benchmark/datavzrd/{species}.{build}.{release}.{datatype}/homer/{macs2_peak_type}/{sample}.render.tsv"
    params:
        extra="",
    wrapper:
        "v3.3.3/utils/datavzrd"
