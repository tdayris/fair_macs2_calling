# -*- coding: utf-8 -*-

"""Snakemake wrapper for MultiQC configuration file"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import yaml

from typing import Any

default_config: dict[str, Any] = {
    "title": "Macs2 Peak-Calling quality control report",
    "subtitle": "Produced on raw fastq recieved from sequencer",
    "intro_text": (
        "This pipeline building this report expects "
        "user-provided informations about sequencing "
        "protocol, samples organism, and wet-lab "
        "experimental design."
    ),
    "report_comment": (
        "This report was generated using: "
        "https://github.com/tdayris/fair_macs2_calling"
    ),
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "custom_logo": snakemake.input[0],
    "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
    "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
    "report_header_info": [
        {"Contact E-mail": "bigr@gustaveroussy.fr"},
        {"Application type": "Short-gapped reads"},
        {"Project Type": "Peak-Calling"},
    ],
    "software_versions": {
        "Quality controls": {
            "fastqc": "1.12.1",
            "fastq_screen": "0.15.3",
            "bowtie2": "1.3.1",
            "multiqc": "1.20.0",
        },
        "Mapping": {
            "bowtie2": "2.5.3",
            "sambamba": "1.0",
            "samtools": "1.19.2",
            "picard": "3.1.1",
            "rseqc": "5.0.3",
            "fastp": "0.23.4",
            "ngsderive": "3.3.2",
            "goleft": "0.2.4",
        },
        "PeakCalling": {
            "macs2": "2.2.9.1",
            "homer": "4.11",
        },
        "Pipeline": {
            "snakemake": "8.5.3",
            "fair_macs2_calling": "2.1.0",
            "fair_bowtie2_mapping": "3.3.2",
            "fair_fastqc_multiqc": "2.2.6",
            "fair_genome_indexer": "3.4.4",
        },
    },
    "disable_version_detection": True,
    "run_modules": [
        "fastqc",
        "fastq_screen",
        "fastp",
        "bowtie2",
        "samtools",
        "picard",
        "rseqc",
        "ngsderive",
        "goleft_indexcov",
        "macs2",
        "homer",
    ],
    "report_section_order": {
        "fastq_screen": {"order": 1000},
        "ngsderive": {"order": 950},
        "fastqc": {"order": 900},
        "fastp": {"order": 890},
        "bowtie2": {"order": 880},
        "picard": {"order": 870},
        "samtools": {"order": 860},
        "rseqc": {"order": 850},
        "goleft_indexcov": {"order": 840},
        "macs2": {"order": 830},
        "homer": {"order": 820},
        "software_versions": {"order": -1000},
    },
}

config: dict[str, Any] | None = snakemake.params.get("extra", None)
if config is None:
    config = default_config.copy()

with open(str(snakemake.output[0]), "w") as out_yaml_stream:
    out_yaml_stream.write(yaml.dump(config, default_flow_style=False))
