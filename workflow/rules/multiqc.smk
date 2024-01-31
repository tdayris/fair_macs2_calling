rule multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC_PeakCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "PeakCalling",
            },
        ),
        "results/QC/MultiQC_PeakCalling_data.zip",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp",
    params:
        extra=config.get("params", {}).get(
            "multiqc", "--zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"
        ),
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.3.3/bio/multiqc"
