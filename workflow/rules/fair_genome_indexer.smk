module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="2.2.0")
    config:
        {"genomes": config["genomes"]}


use rule * from fair_genome_indexer as fair_genome_indexer_*
