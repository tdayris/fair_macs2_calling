module fair_bowtie2_mapping:
    snakefile:
        github(
            "tdayris/fair_bowtie2_mapping",
            path="workflow/Snakefile",
            tag="4.4.6",
        )
    config:
        {
            **config,
            "load_fair_genome_indexer": False,
            "load_fair_fastqc_multiqc": False,
        }


use rule * from fair_bowtie2_mapping
