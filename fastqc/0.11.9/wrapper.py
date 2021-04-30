"""Snakemake wrapper for fastqc."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from os import path
import re
from tempfile import TemporaryDirectory

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
paired = snakemake.params.get("paired", "false")


def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)
    # Remove file extension(s) (similar to the internal fastqc approach)
    base = re.sub("\\.gz$", "", base)
    base = re.sub("\\.bz2$", "", base)
    base = re.sub("\\.txt$", "", base)
    base = re.sub("\\.fastq$", "", base)
    base = re.sub("\\.fq$", "", base)
    base = re.sub("\\.fq.gz$", "", base)
    base = re.sub("\\.fastq$.gz$", "", base)
    base = re.sub("\\.sam$", "", base)
    base = re.sub("\\.bam$", "", base)

    return base


# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
if (paired == "true"):
    with TemporaryDirectory() as tempdir:
        shell(
            "fastqc {snakemake.params} -t {snakemake.threads} "
            "--outdir {tempdir:q} {snakemake.input[0]:q} {snakemake.input[1]:q}"
            " {log}"
        )
        
        file_base_R2 = basename_without_ext(snakemake.input[1])
        output_base_R2 = re.sub("R1", "R2", snakemake.output.html)
        html_path_R2 = path.join(tempdir, file_base_R2 + "_R2_fastqc.html")
        zip_path_R2 = path.join(tempdir, file_base_R2 + "_R2_fastqc.zip")
    
        if snakemake.output.html != html_path_R2:
            shell("mv {html_path:q} {output_base_R2:q}")

        if snakemake.output.zip != zip_path_R2:
            shell("mv {zip_path:q} {output_base_R2.zip:q}")


with TemporaryDirectory() as tempdir:
    shell(
        "fastqc {snakemake.params} -t {snakemake.threads} "
        "--outdir {tempdir:q} {snakemake.input[0]:q}"
        " {log}"
        )
    
        
    # Move outputs into proper position.
    output_base_R1 = basename_without_ext(snakemake.input[0])
    html_path_R1 = path.join(tempdir, output_base_R1 + "_R1_fastqc.html")
    zip_path_R1 = path.join(tempdir, output_base_R1 + "_R1_fastqc.zip")

    if snakemake.output.html != html_path_R1:
        shell("mv {html_path:q} {snakemake.output.html:q}")

    if snakemake.output.zip != zip_path_R1:
        shell("mv {zip_path:q} {snakemake.output.zip:q}")


