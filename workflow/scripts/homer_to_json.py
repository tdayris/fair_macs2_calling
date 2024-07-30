# coding: utf-8

import json
import logging
import pandas

# import snakemake

numeric = int | float


def read_homer_table(path: str) -> pandas.DataFrame:
    """
    Load Homer table in memory, rename columns, define
    index and return a pandas DataFrame

    Parameters:
    path (str): Path to home table

    Return (pandas.DataFrame):
    Loaded table
    """
    logging.debug(f"Reading homer at: {path}")
    df: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=path,
        sep="\t",
        header=0,
    )

    # Consider cases where CpG and GC % are computed
    try:
        df.columns = [
            "PeakID",
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "PeakScore",
            "FocusRatioRegionSize",
            "CompleteAnnotation",
            "DetailedAnnotation",
            "DistanceTSS",
            "NearestPromoterID",
            "EntrezID",
            "NearestUnigene",
            "NearestRefseq",
            "NearestEnsembl",
            "GeneName",
            "GeneAlias",
            "GeneDescription",
            "OriginalGeneType",
            "AvgWig",
        ]
    except ValueError:
        df.columns = [
            "PeakID",
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "PeakScore",
            "FocusRatioRegionSize",
            "CompleteAnnotation",
            "DetailedAnnotation",
            "DistanceTSS",
            "NearestPromoterID",
            "EntrezID",
            "NearestUnigene",
            "NearestRefseq",
            "NearestEnsembl",
            "GeneName",
            "GeneAlias",
            "GeneDescription",
            "OriginalGeneType",
            "PctCpG",
            "pctGC",
            "AvgWig",
        ]
    df.set_index("PeakID", inplace=True)

    # Fix attribute error
    df.dropna(subset=["DistanceTSS"], inplace=True)

    df["Annotations"] = [
        annotation.split(" (")[0].capitalize()
        for annotation in df["CompleteAnnotation"]
    ]
    df["GeneTypes"] = df["OriginalGeneType"].str.lower().replace("rna", "RNA")

    df.sort_values(by=["Chromosome", "Start", "End"], ascending=True)
    logging.debug(df.head())
    return df


def pandas_to_json(df: pandas.DataFrame, sample_id: str) -> dict[str, str | numeric]:
    """
    Summarize Homer annotations of the given sample.

    Parameters:
    df          (pandas.DataFrame): Loaded Homer annotations
    sample_id   (str)             : Sample identifier

    Return (dict[str, str | numeric]):
    Summarized Homer information
    """
    logging.debug("Summurizing:")
    logging.debug(df.head())
    return {
        "name": sample_id,
        "peak_per_chr": dict(
            [key, int(value)] for key, value in df["Chromosome"].value_counts().items()
        ),
        "peak_per_gene_type": dict(
            [key, int(value)] for key, value in df["GeneTypes"].value_counts().items()
        ),
        "peak_per_annotation": dict(
            [key, int(value)] for key, value in df["Annotations"].value_counts().items()
        ),
        "pct_peak_per_chr": dict(
            [key, round(number=value * 100, ndigits=1)]
            for key, value in df["Chromosome"].value_counts(normalize=True).items()
        ),
        "pct_peak_per_gene_type": dict(
            [key, round(number=value * 100, ndigits=1)]
            for key, value in df["GeneTypes"].value_counts(normalize=True).items()
        ),
        "pct_peak_per_annotation": dict(
            [key, round(number=value * 100, ndigits=1)]
            for key, value in df["Annotations"].value_counts(normalize=True).items()
        ),
    }


def save_json(data: dict[str, str | numeric], path: str) -> None:
    """
    Save summary as json data

    Parameters:
    data (dict[str, str | numeric]) : Summarized data
    path (str)                      : Path to output json

    Return: None
    """
    logging.debug(f"Saving summary to {path}")
    with open(path, "w") as json_stream:
        json_obj = json.dumps(data)
        json_stream.write(json_obj)


def main(
    homer_table_path: str, json_output_path: str, sample_id: str, tsv_output_path: str
) -> None:
    """
    Main function: (1) Load homer, (2) Summarize, (3) Save

    Parameters:
    homer_table_path (str): Path to input homer table
    json_output_path (str): Path to output json data
    tsv_output_path  (str): Path to output TSV data

    Return: None
    """
    homer_data: pandas.DataFrame = read_homer_table(path=homer_table_path)
    summary: dict[str, str | numeric] = pandas_to_json(
        df=homer_data, sample_id=sample_id
    )
    save_json(data=summary, path=json_output_path)
    homer_data.to_csv(tsv_output_path, sep="\t", header=True, index=False)


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    try:
        main(
            homer_table_path=snakemake.input[0],
            json_output_path=snakemake.output["json"],
            sample_id=snakemake.wildcards["sample"],
            tsv_output_path=snakemake.output["tsv"],
        )
    except Exception as e:
        logging.error(e)
        raise e
