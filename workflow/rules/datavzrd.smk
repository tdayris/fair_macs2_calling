rule datavzrd_homer_yaml:
    input:
        "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv",
    output:
        temp("tmp/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.yaml"),
    log:
        "logs/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.config.log",
    benchmark:
        "benchmark/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.config.tsv",
    params:
        extra="",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/build_datavzrd_yaml.py"


rule datavzrd_homer_render:
    input:
        config="tmp/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.yaml",
        table="results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv",
    output:
        report(
            directory("results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}_reports/{sample}"),
            htmlindex="index.html",
            caption="../report/datavzrd_homer.rst",
            category="Peak Calling",
            labels={
                "species": "{species}.{build}.{release}",
                "sample": "{sample}",
                "peaks": "{macs2_peak_type}"
            },
        ),
    log:
        "logs/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.render.log",
    benchmark:
        "benchmark/datavzrd/{spacies}.{build}.{release}.{datatype}/homer/{macs2_preak_type}/{sample}.render.tsv",
    params:
        extra="",
    wrapper:
        f"{snakemake_wrappers_version}/utils/datavzrd"