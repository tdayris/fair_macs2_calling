rule deeptools_multibigwig_summary:
    input:
        bed="tmp/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed",
        bw=expand(
            "results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
            sample=samples.sample_id,
            allow_missing=True,
        ),
    output:
        temp(
            "tmp/deeptools/multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz"
        ),
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
    shell:
        "multiBigwigSummary BED-file "
        "--bwfiles {input.bw} "
        "--outFileName {output} "
        "--BED {input.bed} "
        "--numberOfProcessors {threads} "
        "{params.extra} "
        "> {log} 2>&1"


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
    log:
        "logs/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=config.get("params", {}).get("deeptools", {}).get("plot_pca", ""),
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--plotFile {output.png} "
        "--outFileNameData {output.tab} "
        "{params.extra} "
        "> {log} 2>&1"


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
    shell:
        "plotCorrelation "
        "--corData {input} "
        "--plotFile {output.png} "
        "--outFileCorMatrix {output.tab} "
        "{params.extra} "
        "> {log} 2>&1"
