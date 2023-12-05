module fair_bowtie2_mapping:
    snakefile:
        github("tdayris/fair_bowtie2_mapping", path="workflow/Snakefile", tag="2.2.5")
    config:
        {
            "samples": config["samples"],
            "params": config["params"],
            "load_fair_genome_indexer": False,
            "genomes": config["genomes"],
        }


use rule * from fair_bowtie2_mapping as fair_bowtie2_mapping_*


use rule fair_bowtie2_mapping_multiqc_report from fair_bowtie2_mapping with:
    input:
        unpack(get_multiqc_macs2_peakcalling_report_input),
