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
report_dir = Path("report").resolve()
qc_dir = Path("QC").resolve()
mykrobe_dir = Path("mykrobe").resolve()
captions = config["captions"]
reference = config["h37rv"]
pileup_dir = Path("pileup").resolve()
calls_dir = Path("calls").resolve()
filters = config["filters"]
lineage_dir = Path("lineage").resolve()

# ======================================================
# Global functions and variables
# ======================================================
output_files = set()
output_files.add(lineage_dir / "lemur.lineage.csv")
output_files.add(qc_dir / "qc.html")
output_files.add(qc_dir / "lemur.krona.html")
output_files.add(mykrobe_dir / "lemur.dst.json")


# ======================================================
# Rules
# ======================================================
localrules:
    all,
    combine_basecall_fastq,


report: report_dir / "workflow.rst"


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


rule qc_plots:
    input:
        summary=rules.basecall.output.summary,
    output:
        report=report(
            qc_dir / "qc.html", category="Quality Control", caption=captions["pycoqc"],
        ),
    threads: 1
    resources:
        mem_mb=int(GB),
    params:
        options="--report_title 'All basecalled lemur nanopore reads'",
    container:
        containers["pycoqc"]
    log:
        rule_log_dir / "qc_plots.log",
    shell:
        """
        pycoQC {params.options} -f {input.summary} -o {output.report} 2> {log}
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
        alignment=rules.map_reads_to_decontam_db.output.bam,
    output:
        index=mapped_dir / "lemur.all.sorted.bam.bai",
    threads: 4
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["conda"]
    conda:
        envs["aln_tools"]
    log:
        rule_log_dir / "index_mapped_reads.log",
    shell:
        "samtools index -@ {threads} {input.alignment} 2> {log}"


rule filter_contamination:
    input:
        bam=rules.map_reads_to_decontam_db.output.bam,
        index=rules.index_mapped_reads.output.index,
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


rule generate_krona_input:
    input:
        bam=rules.map_reads_to_decontam_db.output.bam,
        index=rules.index_mapped_reads.output.index,
        metadata=rules.filter_contamination.input.metadata,
    output:
        krona_input=qc_dir / "lemur.krona.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(0.5 * GB) * attempt,
    container:
        containers["conda"]
    conda:
        envs["generate_krona_input"]
    params:
        script=scripts["generate_krona_input"],
        extras="--ignore-secondary",
    log:
        rule_log_dir / "generate_krona_input.log",
    shell:
        """
        python {params.script} {params.extras} \
            -i {input.bam} -m {input.metadata} -o {output.krona_input} 2> {log}
        """


rule plot_sample_composition:
    input:
        tsv=rules.generate_krona_input.output.krona_input,
    output:
        chart=report(
            qc_dir / "lemur.krona.html",
            category="Quality Control",
            caption=captions["krona"],
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * GB,
    container:
        containers["krona"]
    log:
        rule_log_dir / "plot_sample_composition.log",
    shell:
        "ktImportText {input.tsv} -o {output.chart} &> {log}"


rule predict_resistance:
    input:
        reads=rules.extract_decontaminated_reads.output.reads,
    output:
        json=report(
            mykrobe_dir / "lemur.dst.json",
            category="AMR Prediction",
            caption=captions["mykrobe"],
        ),
    threads: 1
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["mykrobe"]
    params:
        sample="lemur",
        species="tb",
        options="--force --ont --format json",
    log:
        rule_log_dir / "predict_resistance.log",
    shell:
        """
        mykrobe predict {params.options}  -t {threads} --seq {input.reads} \
            --output {output.json} {params.sample} {params.species} 2> {log}
        """


rule index_reference:
    input:
        reference["genome"],
    output:
        reference["genome"] + ".fai",
    log:
        rule_log_dir / "index_reference.log",
    container:
        containers["conda"]
    wrapper:
        "0.63.0/bio/samtools/faidx"


rule map_reads_to_h37rv:
    input:
        target=rules.index_reference.input[0],
        query=rules.extract_decontaminated_reads.output.reads,
    output:
        bam=mapped_dir / "lemur.filtered.h37rv.sorted.bam",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: int(8 * GB) * attempt,
    container:
        containers["conda"]
    conda:
        envs["aln_tools"]
    params:
        minimap_options=" ".join(
            ["-a", "-L", "--sam-hit-only", "--secondary=no", "-2"]
        ),
        minimap_preset="map-ont",
        sort_options="-O bam",
    log:
        rule_log_dir / "map_reads_to_h37rv.log",
    shell:
        """
        (minimap2 {params.minimap_options} \
            -x {params.minimap_preset} \
            -t {threads} \
            {input.target} {input.query} | \
        samtools sort -@ {threads} -o {output.bam}) 2> {log}
        """


rule pileup:
    input:
        index=rules.index_reference.output[0],
        ref=rules.index_reference.input[0],
        alignments=rules.map_reads_to_h37rv.output.bam,
    output:
        pileup=pileup_dir / "lemur.pileup.bcf",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: int(32 * GB) * attempt,
    params:
        options="--ignore-overlaps -O b -Q 7",
    log:
        rule_log_dir / "pileup.log",
    container:
        containers["conda"]
    wrapper:
        "0.64.0/bio/bcftools/mpileup"


rule call_snps:
    input:
        pileup=rules.pileup.output.pileup,
    output:
        calls=calls_dir / "lemur.snps.bcf",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: int(GB) * attempt,
    params:
        caller="-m",
        options=" ".join(["--ploidy 1", "-O b", "--skip-variants indels"]),
    log:
        rule_log_dir / "call_snps.log",
    container:
        containers["conda"]
    wrapper:
        "0.64.0/bio/bcftools/call"


rule filter_snps:
    input:
        vcf=rules.call_snps.output.calls,
    output:
        vcf=calls_dir / "lemur.snps.filtered.bcf",
    threads: 1
    resources:
        mem_mb=int(GB),
    params:
        options=" ".join(
            [
                "--hist",
                "--verbose",
                "--overwrite",
                f"-d {filters['min_depth']}",
                f"-D {filters['max_depth']}",
                f"-q {filters['min_qual']}",
                f"-s {filters['min_strand_bias']}",
                f"-b {filters['min_bqb']}",
                f"-m {filters['min_mqb']}",
                f"-r {filters['min_rpb']}",
                f"-V {filters['min_vdb']}",
                f"-G {filters['max_sgb']}",
            ]
        ),
        script=scripts["filter_snps"],
    log:
        rule_log_dir / "filter_snps.log",
    container:
        containers["conda"]
    conda:
        envs["filter_snps"]
    shell:
        """
        python {params.script} {params.options} \
            -i {input.vcf} \
            -o {output.vcf} 2> {log}
        """


rule assign_lineage:
    input:
        vcf=rules.filter_snps.output.vcf,
        panel=config["lineage_panel"],
    output:
        assignments=report(
            lineage_dir / "lemur.lineage.csv",
            category="Lineage",
            caption=captions["lineage"],
        ),
    threads: 1
    resources:
        mem_mb=GB,
    container:
        containers["conda"]
    conda:
        envs["assign_lineages"]
    params:
        script=scripts["assign_lineages"],
        default_lineage=config["default_lineage"], # the name given to samples with no hits in the panel
        max_het=1,
        max_alt_lineages=1,
        ref_lineage_position=config["ref_lineage_position"],
        extras="--verbose",
    log:
        rule_log_dir / "assign_lineage.log",
    shell:
        """
        python {params.script} --input {input.vcf} \
            --panel {input.panel} \
            --output {output.assignments} \
            --default-lineage {params.default_lineage} \
            --max-het {params.max_het} \
            --ref-lineage-position {params.ref_lineage_position} \
            --max-alt-lineages {params.max_alt_lineages} {params.extras} 2> {log}
        """
