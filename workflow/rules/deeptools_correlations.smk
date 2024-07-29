rule fair_macs2_calling_deeptools_multibigwig_summary:
    input:
        unpack(get_deeptools_multibigwig_summary_input),
    output:
        npz=temp(
            "tmp/fair_macs2_calling_deeptools_multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_deeptools_multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_multibigwig_summary",
            default="",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_multi_bigwig_summary_wrapper.py"


rule fair_macs2_calling_deeptools_plot_pca:
    input:
        "tmp/fair_macs2_calling_deeptools_multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz",
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
            "tmp/fair_macs2_calling_deeptools_plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_deeptools_plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_plot_pca",
            default="--ntop 1000",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_pca_wrapper.py"


rule fair_macs2_calling_deeptools_plot_correlation:
    input:
        "tmp/fair_macs2_calling_deeptools_multibigwig_summary/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.npz",
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
            "tmp/fair_macs2_calling_deeptools_plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_deeptools_plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_plot_correlation",
            default="--whatToPlot heatmap --corMethod spearman --skipZeros --plotNumbers --colorMap RdYlBu",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_correlation_wrapper.py"


rule fair_macs2_calling_deeptools_plot_enrichment:
    input:
        bams=branch(
            lookup_config(dpath="params/make_sieve", default=False, config=config),
            then=expand(
                "tmp/fair_bowtie2_mapping_sambamba_sort_sieve/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.bam",
                sample=lookup(
                    query="species == '{species}' & release == '{release}' & build == '{build}'",
                    within=samples,
                ),
                datatype="{datatype}",
            ),
            otherwise=expand(
                "results/{sample.species}.{sample.build}.{sample.release}.{datatype}/Mapping/{sample.sample_id}.bam",
                sample=lookup(
                    query="species == '{species}' & release == '{release}' & build == '{build}'",
                    within=samples,
                ),
                datatype="{datatype}",
            ),
        ),
        bais=branch(
            lookup_config(dpath="params/make_sieve", default=False, config=config),
            then=expand(
                "tmp/fair_bowtie2_mapping_sambamba_sort_sieve/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.bam.bai",
                sample=lookup(
                    query="species == '{species}' & release == '{release}' & build == '{build}'",
                    within=samples,
                ),
                datatype="{datatype}",
            ),
            otherwise=expand(
                "results/{sample.species}.{sample.build}.{sample.release}.{datatype}/Mapping/{sample.sample_id}.bam.bai",
                sample=lookup(
                    query="species == '{species}' & release == '{release}' & build == '{build}'",
                    within=samples,
                ),
                datatype="{datatype}",
            ),
        ),
        bed="tmp/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed",
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
            "tmp/fair_macs2_calling_deeptools_plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 45) + (1024 * 20),
        runtime=lambda wildcards, attempt: attempt * 120 + 60,
        tmpdir="tmp",
    log:
        "logs/fair_macs2_calling_deeptools_plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling_deeptools_plot_enrichment/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_deeptools_plot_enrichment",
            default="--ignoreDuplicates --minMappingQuality 30 --samFlagExclude 4 --smartLabels",
        ),
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/deeptools_plot_enrichment_wrapper.py"
