# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Optional parameters
log: str = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra: str = snakemake.params.get("extra", "")

blacklist: str | None = snakemake.input.get("blacklist")
if blacklist:
    extra += f" --blackListFileName {blacklist} "

out_raw_counts: str | None = snakemake.output.get("out_raw_counts")
if out_raw_counts:
    extra += f" --outRawCounts {out_raw_counts} "


shell(
    "multiBigwigSummary BED-file "
    "--bwfiles {snakemake.input.bw} "
    "--outFileName {snakemake.output.npz} "
    "--BED {snakemake.input.bed} "
    "--numberOfProcessors {snakemake.threads} "
    "{extra} {log}"
)
