rule summarize_homer:
    input:
        "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv",
    output:
        json=temp(
            "tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.json"
        ),
        tsv=temp(
            "tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/homer_to_json.py"


rule merge_homer_summaries:
    input:
        unpack(get_merge_homer_summaries_input),
    output:
        annot=temp(
            "tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/Annotations.tsv"
        ),
        chrom=temp(
            "tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/Chromosomes.tsv"
        ),
        genes=temp(
            "tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/GeneTypes.tsv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_homer_summaries.py"


rule plot_homer_summaries_for_regions:
    input:
        summaries="tmp/fair_macs2_calling/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.tsv",
    output:
        catplot=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/{content}.catplot.png",
            caption="../report/catplot.rst",
            category="Peak Calling",
            subcategory="{macs2_peak_type}",
            labels={
                "figure": "catplot",
                "species": "{species}.{build}.{release}",
                "sample": "all",
            },
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 25,
        tmpdir=tmp,
    log:
        "logs/fair_macs2_calling/plot_homer_summaries/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.log",
    benchmark:
        "benchmark/fair_macs2_calling/plot_homer_summaries/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.tsv"
    params:
        content="{content}",
        title=lambda wildcards: f"Percent of {wildcards.macs2_peak_type} falling over Homer {wildcards.content}",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_homer_multisample.py"
