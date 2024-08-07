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


out_raw_counts: str | None = snakemake.output.get("out_raw_counts")
if out_raw_counts:
    extra += f" --outRawCounts {out_raw_counts} "

blacklist: str | None = snakemake.input.get("blacklist")
if blacklist:
    extra += f" --blackListFileName {blacklist} "


shell(
    "plotEnrichment "
    "--bamfiles {snakemake.input.bams} "
    "--BED {snakemake.input.bed} "
    "--plotFile {snakemake.output.png} "
    "--numberOfProcessors {snakemake.threads} "
    "{extra} {log} "
)
