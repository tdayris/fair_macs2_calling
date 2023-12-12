# codinf: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging

# import snakemake
import yaml

from typing import Any, Dict


def write_config(data: Dict[str, Any], path: str) -> None:
    """
    Save dictionary in a yaml-formatted file

    Parameters:
    data    (Dict[str, Any]): datavzrd configuration
    path    (str)           : path to output yaml file

    Return: None
    """
    logging
    with open(path, "w") as yaml_stream:
        yaml.dump(data=data, stream=yaml_stream, default_flow_style=False)


def fill_dict(
    sample_name: str, table_path: str, table_description: str
) -> Dict[str, Any]:
    """
    Build datavzrd dictionary

    Parameters:
    sample_name         (str): Sample identifier
    table_path          (str): Path to the data table
    table_description   (str): Header text in the HTML index table

    Return (Dict[str, Any]):
    Datavzrd configuration
    """
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

if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    try:
        datavzrd_configuration: Dict[str, Any] = fill_dict(
            sample_name=snakemake.wildcards.sample,
            table_path=snakemake.input.summary,
            table_description=table_description,
        )
        write_config(data=datavzrd_configuration, path=snakemake.output.yaml)
    except Exception as e:
        logging.error(e)
        raise e
