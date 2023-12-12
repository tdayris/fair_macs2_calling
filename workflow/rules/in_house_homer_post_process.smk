rule summarize_homer:
    input:
        "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv",
    output:
        json=temp(
            "tmp/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.json"
        ),
        tsv=temp(
            "tmp/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
        ),
    log:
        "logs/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/homer_to_json.py"


rule merge_homer_summaries:
    input:
        unpack(get_merge_homer_summaries_input),
    output:
        annot=temp(
            "tmp/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/Annotations.tsv"
        ),
        chrom=temp(
            "tmp/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/Chromosomes.tsv"
        ),
        genes=temp(
            "tmp/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/GeneTypes.tsv"
        ),
    log:
        "logs/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_homer_summaries.py"


rule plot_homer_summaries_for_regions:
    input:
        summaries="tmp/summarize_homer/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.tsv",
    output:
        catplot=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/{content}.catplot.png",
            caption="../report/catplot.rst",
            category="Peak Calling",
            subcategory="{macs2_peak_type}",
            labels={
                "figure": "catplot",
                "species": "{species}.{build}.{release}",
            },
        ),
        dotplot=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/{content}.dotplot.png",
            caption="../report/dotplot.rst",
            category="Peak Calling",
            subcategory="{macs2_peak_type}",
            labels={
                "figure": "dotplot",
                "species": "{species}.{build}.{release}",
            },
        ),
    log:
        "logs/plot_homer_summaries/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.log",
    benchmark:
        "benchmark/plot_homer_summaries/{species}.{build}.{release}.{datatype}.{macs2_peak_type}/{content}.tsv"
    params:
        content="{content}",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_homer_multisample.py"
