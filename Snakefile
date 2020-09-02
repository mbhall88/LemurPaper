import os
from pathlib import Path
from typing import Dict, Union, List

GB = 1_024
PathLike = Union[str, Path, os.PathLike]


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


containers: Dict[str, PathLike] = config["containers"]
# envs: Dict[str, PathLike] = config["envs"]
# scripts: Dict[str, PathLike] = config["scripts"]
rule_log_dir = Path("logs/stderr").resolve()
fast5_dir = Path(config["fast5_dir"]).resolve()
basecall_dir = Path("basecall").resolve()
guppy_outdir = basecall_dir / "guppy"

# ======================================================
# Global functions and variables
# ======================================================
output_files = set()
output_files.add(basecall_dir / "lemur.all.fq")


# ======================================================
# Rules
# ======================================================
localrules:
    all,
    combine_basecall_fastq,


rule all:
    input:
        output_files,


rule basecall:
    input:
        fast5=fast5_dir,
    output:
        summary=guppy_outdir / "sequencing_summary.txt",
    threads: 2
    resources:
        mem_mb=int(4 * GB),
    container:
        containers["guppy"]
    params:
        num_callers=8,
        save_path=lambda wildcards, output: Path(output.summary).parent,
        device="cuda:all:100%",
        config=config["model_config"],
        extras="--recursive",
    log:
        rule_log_dir / "basecall.log",
    shell:
        """
        guppy_basecaller {params.extras} \
            --input_path {input.fast5} \
            --save_path {params.save_path} \
            --config {params.config} \
            --device {params.device} \
            --num_callers {params.num_callers}  &> {log}
        """


rule combine_basecall_fastq:
    input:
        summary=rules.basecall.output.summary,
    output:
        reads=basecall_dir / "lemur.all.fq",
    threads: 1
    resources:
        mem_mb=int(0.5 * GB),
    params:
        indir=lambda wildcards, input: Path(input.summary).parent,
    log:
        rule_log_dir / "combine_basecall_fastq.log",
    shell:
        """
        awk 1 {params.indir}/*.fastq > {output.reads} 2> {log}
        """
