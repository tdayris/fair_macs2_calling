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


out_file_name_data: str | None = snakemake.input.get("tab")
if out_file_name_data:
    extra += f" --outFileNameData {out_file_name_data} "
    
shell(
    "plotPCA "
    "--corData {input} "
    "--plotFile {output.png} "
    "{extra} {log} "
)