rule multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC.PeakCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="PeakCalling",
            labels={
                "report": "html",
            },
        ),
        "results/QC/MultiQC.PeakCalling_data.zip",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.2.0/bio/multiqc"
