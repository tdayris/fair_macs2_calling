rule homer_annotate_peaks:
    input:
        unpack(get_homer_annotate_peaks_input),
    output:
        annotations=protected(
            ensure(
                "results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv",
                non_empty=True,
            ),
        ),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 6),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    shadow:
        "minimal"
    log:
        "logs/fair_macs2_calling/homer/annotatepeaks/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.log",
    benchmark:
        "benchmark/fair_macs2_calling/homer/annotatepeaks/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.tsv"
    params:
        extra=lambda wildcards: get_homer_annotate_peaks_params(
            wildcards, samples, genomes, config
        )["extra"],
    conda:
        "../envs/homer.yaml"
    script:
        "../scripts/homer_wrapper.py"
