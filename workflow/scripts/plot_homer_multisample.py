# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import csv
import logging
import pandas
import seaborn

# import snakemake

import matplotlib.pyplot as plt


def read_table(path: str) -> pandas.DataFrame:
    """
    Load a CSV/TSV table in memory

    Parameters:
    path (str): Path to the table

    Return (pandas.DataFrame):
    Loaded table
    """
    logging.debug(f"Loading {path=}")
    with open(path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.read(1024))
        table_stream.seek(0)

    return pandas.read_csv(
        filepath_or_buffer=path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        # dtype=str,
    )


def catplot(df: pandas.DataFrame, x: str, out_png: str, title: str) -> None:
    """
    Save a catplot on the disk, using x/hue to highlight categories

    Parameters:
    df      (pandas.DataFrame)  : Summary informations
    x       (str)               : Name of the column of interest
    out_png (str)               : Path to output PNG file
    title   (str)               : Graph title

    Return: None
    """
    logging.debug(f"Saving catplot to {out_png=}")
    seaborn.set_theme(style="whitegrid")

    graph = seaborn.catplot(
        data=df,
        x=x,
        y="Percent",
        hue="Sample_id",
        errorbar="sd",
        palette="tab10",
        kind="bar",
        legend_out=False,
    )

    graph.despine(left=True)
    graph.set_axis_labels("", "Percent of peaks")
    graph.set(title=title)

    graph.set_xticklabels(rotation=90)

    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.cla()
    plt.clf()
    plt.close()


def dotplots(df: pandas.DataFrame, x: str, out_png: str, title: str) -> None:
    """
    Save categorical dotplots on disk

    Parameters:
    df      (pandas.DataFrame)  : Summary informations
    x       (str)               : Name of the column of interest
    out_png (str)               : Path to output PNG file
    title   (str)               : Graph title

    Return: None
    """
    logging.debug(f"Saving dotplots to {out_png=}")
    seaborn.set_theme(style="whitegrid")

    data: pandas.DataFrame = (
        df.pivot(index="Sample_id", columns=x, values="Percent")
        .fillna(0.0)
        .sort_index()
        .round(1)
        .astype(float)
    )
    logging.debug(data.head())

    graph = seaborn.PairGrid(
        data=df,
        x_vars=list(df.columns),
        y_vars=list(df.index),
        height=10,
        aspect=0.25,
    )

    graph.map(
        seaborn.stripplot,
        size=10,
        orient="h",
        jitter=False,
        palette="flate_r",
        linewidth=1,
        edgecolor="w",
    )

    graph.set(xlim=(0, 100), xlabel="Percent", ylabel="")

    titles: list[str] = list(df.columns)
    for ax, title in zip(graph.axes, titles):
        ax.set(title=title)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

    seaborn.despine(left=True, bottom=True)

    plt.tight_layout()
    plt.savefig(out_png, dpi=100)
    plt.cla()
    plt.clf()
    plt.close()


def main(
    summaries_path: str, content: str, out_catplot: str, out_dotplot: str, title: str
) -> None:
    """
    Main function: (1) Load table, (2) Save CatPlot, (3) Save DotPlot

    Parameters:
    summaries_path  (str): Path to the summary table
    content         (str): Column of interest
    out_catplot     (str): Path to catplot PNG file
    out_catplot     (str): Path to dotplot PNG file
    catplot_title   (str): Title of the catplot graph

    Return: None
    """
    summary: pandas.DataFrame = read_table(path=summaries_path)
    logging.debug(summary.head())
    catplot(df=summary, x=content, out_png=out_catplot, title=title)
    if out_dotplot:
        dotplots(df=summary, x=content, out_png=out_dotplot, title=title)


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    try:
        main(
            summaries_path=snakemake.input.summaries,
            content=snakemake.params.content,
            out_catplot=snakemake.output.get("catplot"),
            out_dotplot=snakemake.output.get("dotplot"),
            title=snakemake.params.title,
        )
    except Exception as e:
        logging.error(e)
        raise e
