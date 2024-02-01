rule deeptools_multibigwig_summary:
    input:
        unpack(get_deeptools_multibigwig_summary_input),
    output:
        npz=temp(
            "tmp/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {})
        .get("deeptools", {})
        .get("multibigwig_summary", ""),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_multi_bigwig_summary_wrapper.py"


rule deeptools_plot_pca:
    input:
        "tmp/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz",
    output:
        png=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/PCA.png",
            caption="../report/deeptools_plotpca.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_pca",
                "species": "{species}.{build}.{release}",
            },
        ),
        tab=temp(
            "tmp/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {}).get("deeptools", {}).get("plot_pca", ""),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_pca_wrapper.py"


rule deeptools_plot_correlation:
    input:
        "tmp/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz",
    output:
        png=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/Heatmap.png",
            caption="../report/deeptools_plotcorrelation.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_heatmap",
                "species": "{species}.{build}.{release}",
            },
        ),
        tab=temp(
            "tmp/deeptools/plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/deeptools/plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/deeptools/plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {})
        .get("deeptools", {})
        .get("plot_correlation", "--corMethod spearman --whatToPlot heatmap"),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_correlation_wrapper.py"


rule deeptools_plot_enrichment:
    input:
        unpack(get_deeptools_plotcoverage_input),
    output:
        png=report(
            "results/{species}.{build}.{release}.{datatype}/Graphs/{macs2_peak_type}/Enrichment.png",
            caption="../report/deeptools_plotenrichment.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_enrichment",
                "species": "{species}.{build}.{release}",
            },
        ),
        out_raw_counts=temp(
            "tmp/deeptools/plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/deeptools/plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/deeptools/plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {}).get("deeptools", {}).get("plot_enrichment", ""),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_enrichment_wrapper.py"
