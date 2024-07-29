# -*- coding: utf-8 -*-

"""Snakemake wrapper for MultiQC configuration file"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import yaml
import pandas

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
            "deeptools": "3.5.5",
            "macs2": "2.2.9.1",
            "homer": "4.11",
        },
        "Pipeline": {
            "snakemake": "8.16.0",
            "fair_macs2_calling": "3.0.2",
            "fair_bowtie2_mapping": "4.2.0",
            "fair_fastqc_multiqc": "2.3.5",
            "fair_genome_indexer": "3.8.1",
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
        "deeptools",
        "custom_content",
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
        "deeptools": {"order": 830},
        "macs2": {"order": 820},
        "software_versions": {"order": -1000},
    },
   "custom_data": {
       "BiGR": {
           "id": "bigr_homer_annotation",
           "section_anchor": "bigr_homer_annotation", 
           "section_name": "Homer Annotation",
           "section_href": "http://homer.ucsd.edu/homer/ngs/annotation.html",
           "description": "This graphs describes where annotation falls over genome annotations.",
           "plot_type": "bargraph",
           "pconfig": {
               "id": "barplot_config_only",
               "title": "Number of peaks overlapping genomic annotations",
               "ylab": "Number of peaks",
           },
           "data": {}
       }
   }
}

config: dict[str, Any] | None = snakemake.params.get("extra", None)
if config is None:
    config = default_config.copy()

homer_df = pandas.read_csv(
   snakemake.input["homer_annotations"],
   sep="\t",
   header=0,
   index_col=0,
)
homer_dict: dict[str, float] = {}
for data in homer_df.itertuples():
    if data.Sample_id in homer_dict.keys():
        homer_dict[data.Sample_id][data.Index] = data.Peaks
    else:
        homer_dict[data.Sample_id] = {data.Index: data.Peaks}
    
config["custom_data"]["BiGR"]["data"] = homer_dict


with open(str(snakemake.output[0]), "w") as out_yaml_stream:
    out_yaml_stream.write(yaml.dump(config, default_flow_style=False))
