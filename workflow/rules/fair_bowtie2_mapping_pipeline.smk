module fair_bowtie2_mapping:
    snakefile:
        github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="2.3.0")
    config:
        {
            "samples": config.get("samples", "config/samples.csv"),
            "params": {
                "fastp": config.get("params", {}).get("fastp", {}),
                "fastqc": config.get("params", {}).get("fastqc", ""),
                "bowtie2": config.get("params", {}).get("bowtie2", {}),
                "sambamba": {
                    "view": config.get("params", {})
                    .get("sambamba", {})
                    .get("view", "--format 'bam'"),
                    "markdup": config.get("params", {})
                    .get("sambamba", {})
                    .get("markdup", "--overflow-list-size=500000"),
                },
                "picard": config.get("params", {}).get("picard", {}),
                "samtools": config.get("params", {}).get("samtools", ""),
                "multiqc": config.get("params", {}).get("multiqc", "--zip-data-dir"),
            },
            "genomes": config.get("genomes", "config/genomes.csv"),
            "load_fair_genome_indexer": False,
            "load_fair_fastqc_multiqc": False,
        }


use rule * from fair_bowtie2_mapping
