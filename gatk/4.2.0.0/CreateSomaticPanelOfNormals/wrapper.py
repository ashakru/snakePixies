__author__ = "Joanna Krupka"
__copyright__ = "Copyright 2021, Joanna Krupka"
__email__ = "joannaakrupka@gmail.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

gvcfs = list(map("-V {}".format, snakemake.input))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk CreateSomaticPanelOfNormals {extra} "
    "{gvcfs} "
    "-R {snakemake.params.index} "
    "--germline-resource {snakemake.params.gnomad} "
    "-O {snakemake.output} {log}"
)
