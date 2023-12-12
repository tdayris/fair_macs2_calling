# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import json
import logging
import pandas

# import snakemake

numeric = int | float


def read_json_as_dict(path: str) -> dict[str, numeric]:
    """
    Load a json-formatted homer summary into a dictionary

    Parameters:
    path (str): Path to json-formatted homer sommary

    Return (dict[str, numeric]):
    Homer summary
    """
    logging.debug(f"Loading {path}")
    with open(path, "r") as json_stream:
        return json.load(json_stream)


def dict_to_pandas(summaries: list[dict[str, numeric]]) -> dict[str, pandas.DataFrame]:
    """
    Concatenate all summaries in a single pandas DataFrame

    Parameters:
    summaries (list[dict[str, numeric]]): All summaries to concatenate

    Return (dict[str, pandas.DataFrame]):
    The dictionnary of complete summary tables: annotation, gene, and chromosomes
    """
    logging.debug("Merging all summaries")
    merged: dict[str, pandas.DataFrame] = {
        "gene_type": None,
        "annotation": None,
        "chr": None,
    }
    for content in ["gene_type", "annotation", "chr"]:
        logging.debug(f"Working on {content=}")
        for summary in summaries:
            # cast summary as a dataframe
            summary_df: pandas.DataFrame = pandas.DataFrame.from_records(
                [
                    summary[f"pct_peak_per_{content}"],
                    summary[f"peak_per_{content}"],
                ]
            ).transpose()
            if len(summary_df) == 0:
                continue
            summary_df["name"] = summary["name"]
            summary_df.columns = ["Percent", "Peaks", "Sample_id"]
            summary_df.index.name = (
                "Chromosome"
                if content == "chr"
                else ("Annotation" if content == "annotation" else "Gene_Type")
            )
            summary_df.reset_index(inplace=True)

            # Concat dataframe
            merged[content] = pandas.concat(
                [merged[content], summary_df], ignore_index=True
            )
            logging.debug(merged[content].head())

    logging.debug("Finished concatenation:")
    logging.debug(merged)
    return merged


def pandas_to_file(df: pandas.DataFrame, path: str) -> None:
    """
    Save on disk according to the file extension
    """
    logging.debug(f"Saving table to {path=}")
    logging.debug(df.head())
    sep: str = ","
    if path.endswith("tsv"):
        sep = "\t"

    with open(path, "w") as pandas_stream:
        df.to_csv(
            path_or_buf=pandas_stream,
            sep=sep,
            header=True,
            index=False,
        )


def main(
    summaries_paths: list[str],
    out_table_chrom: str,
    out_table_annot: str,
    out_table_genes: str,
) -> None:
    """
    Main function: (1) Load json, (2) build table, (3) save table

    Parameters:
    summaries_paths (list[str]) : Path to every single Homer summaries
    out_table       (str)       : Output table path

    Return: None
    """
    summaries: list[dict[str, numeric]] = [
        read_json_as_dict(path=summary_path) for summary_path in summaries_paths
    ]
    summaries_tables: pandas.DataFrame = dict_to_pandas(summaries=summaries)

    pandas_to_file(df=summaries_tables["chr"], path=out_table_chrom)
    pandas_to_file(df=summaries_tables["annotation"], path=out_table_annot)
    pandas_to_file(df=summaries_tables["gene_type"], path=out_table_genes)


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    try:
        main(
            out_table_chrom=snakemake.output["chrom"],
            out_table_annot=snakemake.output["annot"],
            out_table_genes=snakemake.output["genes"],
            summaries_paths=snakemake.input["summaries"],
        )
    except Exception as e:
        logging.error(e)
        raise e
