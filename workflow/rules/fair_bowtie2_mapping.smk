module fair_bowtie2_mapping:
    snakefile:
        github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="2.2.1")
    config:
        config


use rule * from fair_bowtie2_mapping

use rule multiqc_report with:
    input:
        unpack(get_multiqc_report_peak_input)