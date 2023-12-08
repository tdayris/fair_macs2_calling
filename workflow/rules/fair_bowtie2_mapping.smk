module fair_bowtie2_mapping:
    snakefile:
        github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="2.2.6")
    config:
        {
            "samples": config["samples"],
            "params": config["params"],
            "load_fair_genome_indexer": False,
            "genomes": config["genomes"],
        }


use rule * from fair_bowtie2_mapping as fair_bowtie2_mapping_*
