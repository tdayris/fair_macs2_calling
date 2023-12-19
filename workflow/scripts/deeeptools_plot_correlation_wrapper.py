# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Optional parameters
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra: str = snakemake.params.get("extra", "")

accepted_plot_file_format: list[str] = ["png", "pdf", "svg", "eps", "plotly"]
plot_file_format: str = str(snakemake.output.png).split(".")[-1].lower()
if plot_file_format in accepted_plot_file_format:
    extra += f" --plotFileFormat {plot_file_format} "
else:
    raise ValueError(f"Plot file format sould be one of {accepted_plot_file_format}")


corr_matrix: str | None = snakemake.output.get("tab")
if corr_matrix:
    extra += f" --outFileCorMatrix {corr_matrix} "

method: str | None = snakemake.wildcards.get("method")
if method:
    extra += f" --corMethod {method} "

plot: str | None = snakemake.wildcards.get("plot")
if plot:
    extra += f" --whatToPlot {plot} "


shell(
    "plotCorrelation "
    "--corData {input} "
    "--plotFile {output.png} "
    "{extra} {log}"
)