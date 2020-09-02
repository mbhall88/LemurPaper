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
envs: Dict[str, PathLike] = config["envs"]
scripts: Dict[str, PathLike] = config["scripts"]
rule_log_dir = Path("logs/stderr").resolve()
fast5_dir = Path(config["fast5_dir"]).resolve()
basecall_dir = Path("basecall").resolve()
guppy_outdir = basecall_dir / "guppy"
decontam_dir = Path(config["decontam_db"]).resolve()
mapped_dir = Path("mapped").resolve()
filtered_dir = Path("filtered").resolve()

# ======================================================
# Global functions and variables
# ======================================================
output_files = set()
output_files.add(filtered_dir / "lemur.filtered.fq")


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


rule map_reads_to_decontam_db:
    input:
        index=decontam_dir / "remove_contam.fa.gz.map-ont.mmi",
        query=rules.combine_basecall_fastq.output.reads,
    output:
        bam=mapped_dir / "lemur.all.sorted.bam",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
    params:
        map_options="-aL2 -x map-ont",
        sort_options="-O bam",
    container:
        containers["conda"]
    conda:
        envs["aln_tools"]
    log:
        rule_log_dir / "map_reads_to_decontam_db.log",
    shell:
        """
        (minimap2 {params.map_options} -t {threads} {input.index} {input.query} | \
            samtools sort -@ {threads} -o {output.bam}) 2> {log}
        """


rule index_mapped_reads:
    input:
        rules.map_reads_to_decontam_db.output.bam,
    output:
        mapped_dir / "lemur.all.sorted.bam.bai",
    threads: 4
    params:
        lambda wildcards, threads: f"-@ {threads}",
    container:
        containers["conda"]
    wrapper:
        "0.65.0/bio/samtools/index"


rule filter_contamination:
    input:
        bam=rules.map_reads_to_decontam_db.output.bam,
        index=rules.index_mapped_reads.output[0],
        metadata=decontam_dir / "remove_contam.tsv",
    output:
        keep_ids=filtered_dir / "keep.reads",
        contam_ids=filtered_dir / "contaminant.reads",
        unmapped_ids=filtered_dir / "unmapped.reads",
    threads: 1
    resources:
        mem_mb=GB,
    container:
        containers["conda"]
    conda:
        envs["filter_reads"]
    params:
        script=scripts["filter_reads"],
        extra="--verbose --ignore-secondary",
        outdir=lambda wildcards, output: Path(output.keep_ids).parent,
    log:
        rule_log_dir / "filter_contamination.log",
    shell:
        """
        python {params.script} {params.extra} \
            -i {input.bam} \
            -m {input.metadata} \
            -o {params.outdir} 2> {log}
        """


rule extract_decontaminated_reads:
    input:
        reads=rules.map_reads_to_decontam_db.input.query,
        read_ids=rules.filter_contamination.output.keep_ids,
    output:
        reads=filtered_dir / "lemur.filtered.fq",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    container:
        containers["conda"]
    conda:
        envs["rg"]
    log:
        rule_log_dir / "extract_decontaminated_reads.log",
    shell:
        """
        (paste - - - - < {input.reads} | \
        rg -f {input.read_ids} | tr "\t" "\n") > {output.reads} 2> {log}
        """
