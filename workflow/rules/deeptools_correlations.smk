rule deeptools_multibigwig_summary:
    input:
        unpack(get_deeptools_multibigwig_summary_input),
    output:
        npz=temp(
            "tmp/fair_macs2_calling/deeptools/multibigwig_summary/{species}.{build}.{release}.dna/{macs2_peak_type}.npz"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling/deeptools/multibigwig_summary/{species}.{build}.{release}.dna/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/multibigwig_summary/{species}.{build}.{release}.dna/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_macs2_calling/deeptools/multibigwig_summary",
            within=config,
            default="",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_multi_bigwig_summary_wrapper.py"


rule deeptools_plot_pca:
    input:
        "tmp/fair_macs2_calling/deeptools/multibigwig_summary/{species}.{build}.{release}.dna/{macs2_peak_type}.npz",
    output:
        png=report(
            "results/{species}.{build}.{release}.dna/Graphs/{macs2_peak_type}/PCA.png",
            caption="../report/deeptools_plotpca.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_pca",
                "species": "{species}.{build}.{release}",
            },
        ),
        tab=temp(
            "tmp/fair_macs2_calling/deeptools/plot_pca/{species}.{build}.{release}.dna/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling/deeptools/plot_pca/{species}.{build}.{release}.dna/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/plot_pca/{species}.{build}.{release}.dna/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_macs2_calling/deeptools/plot_pca",
            within=config,
            default="--ntop 1000",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_pca_wrapper.py"


rule deeptools_plot_correlation:
    input:
        "tmp/fair_macs2_calling/deeptools/multibigwig_summary/{species}.{build}.{release}.dna/{macs2_peak_type}.npz",
    output:
        png=report(
            "results/{species}.{build}.{release}.dna/Graphs/{macs2_peak_type}/Heatmap.png",
            caption="../report/deeptools_plotcorrelation.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_heatmap",
                "species": "{species}.{build}.{release}",
            },
        ),
        tab=temp(
            "tmp/fair_macs2_calling/deeptools/plot_correlation/{species}.{build}.{release}.dna/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling/deeptools/plot_correlation/{species}.{build}.{release}.dna/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/plot_correlation/{species}.{build}.{release}.dna/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_macs2_calling/deeptools/plot_correlation",
            within=config,
            default="--whatToPlot heatmap --corMethod spearman --skipZeros --plotNumbers --colorMap RdYlBu",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_correlation_wrapper.py"


rule deeptools_plot_enrichment:
    input:
        unpack(get_deeptools_plotcoverage_input),
    output:
        png=report(
            "results/{species}.{build}.{release}.dna/Graphs/{macs2_peak_type}/Enrichment.png",
            caption="../report/deeptools_plotenrichment.rst",
            subcategory="{macs2_peak_type}",
            category="Correlation",
            labels={
                "figure": "plot_enrichment",
                "species": "{species}.{build}.{release}",
            },
        ),
        out_raw_counts=temp(
            "tmp/fair_macs2_calling/deeptools/plot_enrichment/{species}.{build}.{release}.dna/{macs2_peak_type}.tab"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling/deeptools/plot_enrichment/{species}.{build}.{release}.dna/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/deeptools/plot_enrichment/{species}.{build}.{release}.dna/{macs2_peak_type}.tsv"
    params:
        extra=dlookup(
            dpath="params/fair_macs2_calling/deeptools/plot_enrichment",
            within=config,
            default="--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --smartLabels",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_enrichment_wrapper.py"
