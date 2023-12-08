# codinf: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import snakemake
import yaml

from typing import Any, Dict

def write_config(data: Dict[str, Any], path: str) -> None:
    with open(path, "w") as yaml_stream:
        yaml.dump(data=data, stream=yaml_stream, default_flow_style=False)


def fill_dict(sample_name: str, table_path: str, table_description: str) -> Dict[str, Any]:
    return {
        "name": f"Annotations of peaks detected by Macs2, and annotated by Homer, over the sample {sample_name}",
        "datasets": {
            "homer": {
                "path": table_path,
                "separator": "\t",
                "headers": 1,
                "offer-excel": False,
            },
        },
        "views": {
            "homer": {
                "dataset": "homer",
                "desc": table_description,
            },
        },
    }


table_description: str = f"""
# Homer annotation

Please find below the Homer annotations for the sample {snakemake.wildcards.sample}.
"""