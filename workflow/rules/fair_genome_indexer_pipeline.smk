module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="2.3.2")
    config:
        {"genomes": config.get("genomes", "config/genomes.csv")}


use rule * from fair_genome_indexer
