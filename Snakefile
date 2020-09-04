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
resource_dir = Path("resources").resolve()
consensus_dir = Path("consensus").resolve()
distance_dir = Path("distance").resolve()
phylo_dir = Path("phylo").resolve()
other_consensuses: List[PathLike] = []
other_consensus_dir = Path(config["other_consensus_dir"]).resolve()
for sample in config["other_consensuses"]:
    other_consensuses.append(other_consensus_dir / f"{sample}.consensus.fa")

# ======================================================
# Global functions and variables
# ======================================================
output_files = set()
output_files.add(lineage_dir / "lemur.lineage.csv")
output_files.add(qc_dir / "qc.html")
output_files.add(qc_dir / "lemur.krona.html")
output_files.add(mykrobe_dir / "lemur.dst.json")
output_files.add(distance_dir / "heatmap.html")
output_files.add(phylo_dir / "lemur.tree")
output_files.add(distance_dir / "non_lemur.dotplot.html")
output_files.add(distance_dir / "non_lemur.heatmap.nanopore.html")
output_files.add(distance_dir / "non_lemur.heatmap.illumina.html")


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


rule download_lemur_assembly:
    output:
        asm=resource_dir / "GCA_004024665.1.fa.gz",
    resources:
        mem_mb=int(0.3 * GB),
    params:
        url=config["lemur_assembly_url"],
    log:
        rule_log_dir / "download_lemur_assembly.log",
    shell:
        "wget {params.url} -O {output.asm} 2> {log}"


rule add_lemur_to_decontam_db:
    input:
        asm=rules.download_lemur_assembly.output.asm,
    output:
        db=resource_dir / "remove_contam_with_lemur.fa.gz",
        metadata=resource_dir / "remove_contam_with_lemur.tsv",
    threads: 1
    resources:
        mem_mb=int(GB),
    log:
        rule_log_dir / "add_lemur_to_decontam_db.log",
    params:
        metadata=decontam_dir / "remove_contam.tsv",
        db=decontam_dir / "remove_contam.fa.gz",
        name="Lemur",
        contamination_code="1",
        seqid_pattern=r"'^>(?P<id>[\w\.]+)\s.*$'",
        replace_with="'$id'",
        rg_opts="-uu -o -N -z",
    container:
        containers["conda"]
    conda:
        envs["rg"]
    shell:
        """
        cat {params.db} {input.asm} > {output.db} 2> {log}
        (rg {params.rg_opts} -r {params.replace_with} {params.seqid_pattern} {input.asm} | \
          awk '{{ print "{params.name}\t{params.contamination_code}\t"$1 }}' | \
          cat {params.metadata} -) > {output.metadata} 2>> {log}
        """


rule map_reads_to_decontam_db:
    input:
        target=rules.add_lemur_to_decontam_db.output.db,
        query=rules.combine_basecall_fastq.output.reads,
    output:
        bam=mapped_dir / "lemur.all.sorted.bam",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(40 * GB),
    params:
        map_options="-aL2 -x map-ont -I 8G",
        sort_options="-O bam",
    container:
        containers["conda"]
    conda:
        envs["aln_tools"]
    log:
        rule_log_dir / "map_reads_to_decontam_db.log",
    shell:
        """
        (minimap2 {params.map_options} -t {threads} {input.target} {input.query} | \
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
        metadata=rules.add_lemur_to_decontam_db.output.metadata,
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


rule mask_h37rv:
    input:
        genome=reference["genome"],
        mask=reference["mask"],
    output:
        masked_genome=resource_dir / "h37rv.masked.fa",
    threads: 1
    resources:
        mem_mb=int(GB),
    container:
        containers["bedtools"]
    log:
        rule_log_dir / "mask_h37rv.log",
    shell:
        """
        bedtools maskfasta -fi {input.genome} -bed {input.mask} \
            -fo {output.masked_genome} 2> {log}
        """


rule generate_consensus:
    input:
        mask=rules.mask_h37rv.input.mask,
        ref_fasta=reference["genome"],
        vcf=rules.filter_snps.output.vcf,
    output:
        fasta=consensus_dir / "lemur.consensus.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(GB) * attempt,
    container:
        containers["conda"]
    conda:
        envs["consensus"]
    params:
        options=" ".join(
            ["--verbose", "--ignore all", "--sample-id lemur", "--het-default none"]
        ),
        script=scripts["consensus"],
    log:
        rule_log_dir / "generate_consensus.log",
    shell:
        """
        python {params.script} {params.options} \
            -i {input.vcf} \
            -f {input.ref_fasta} \
            -m {input.mask} \
            -o {output.fasta} 2> {log}
        """


def infer_non_lemur_consensus_path(wildcards) -> List[str]:
    paths = []
    for sample in config["other_consensuses"]:
        caller = "bcftools" if wildcards.tech == "nanopore" else "compass"
        path = other_consensus_dir / f"{caller}/madagascar/{sample}.consensus.fa"
        paths.append(str(path))
    return paths


rule aggregate_non_lemur_consensus:
    input:
        non_lemur_consensuses=infer_non_lemur_consensus_path,
    output:
        fasta=consensus_dir / "non_lemur.consensus.{tech}.fa",
    threads: 1
    resources:
        mem_mb=int(GB),
    log:
        rule_log_dir / "aggregate_non_lemur_consensus/{tech}.log",
    shell:
        "awk 1 {input.non_lemur_consensuses} > {output.fasta} 2> {log}"


rule non_lemur_snp_distance:
    input:
        fasta=rules.aggregate_non_lemur_consensus.output.fasta,
    output:
        matrix=distance_dir / "non_lemur.matrix.{tech}.csv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(GB) * attempt,
    container:
        containers["snp-dists"]
    log:
        rule_log_dir / "non_lemur_snp_distance/{tech}.log",
    params:
        options="-c",
    shell:
        """
        snp-dists {params.options} {input.fasta} 2> {log} > {output.matrix}
        """


rule plot_non_lemur_distance_matrix:
    input:
        matrix=rules.non_lemur_snp_distance.output.matrix,
    output:
        plot=report(
            distance_dir / "non_lemur.heatmap.{tech}.html",
            category="Distance",
            caption=captions["distance_matrix"],
        ),
    params:
        script=scripts["plot_distance_matrix"],
        options=" ".join(
            ["--delim ,", "--title 'Non Lemur Samples {tech} Pairwise distance'",]
        ),
    threads: 1
    resources:
        mem_mb=int(GB),
    container:
        containers["conda"]
    conda:
        envs["plot_distance_matrix"]
    shell:
        """
        python {params.script} {params.options} -i {input.matrix} -o {output.plot}
        """


rule dotplot_non_lemur_samples:
    input:
        x_matrix=distance_dir / "non_lemur.matrix.illumina.csv",
        y_matrix=distance_dir / "non_lemur.matrix.nanopore.csv",
    output:
        plot=report(
            distance_dir / "non_lemur.dotplot.html",
            caption=report_dir / "dotplot.rst",
            category="Distance",
        ),
    threads: 1
    resources:
        mem_mb=int(2 * GB),
    container:
        containers["conda"]
    conda:
        envs["dotplot"]
    params:
        script=scripts["dotplot"],
        xname="illumina",
        yname="nanopore",
        options=" ".join(
            [
                "--title 'Pairwise SNP distances for illumina and nanopore calls'",
                "--delim ,",
            ]
        ),
    log:
        rule_log_dir / "dotplot_non_lemur_samples.log",
    shell:
        """
        python {params.script} {params.options} \
            -X {params.xname} -Y {params.yname} \
            -x {input.x_matrix} -y {input.y_matrix} \
            -o {output.plot} 2> {log}
        """


rule aggregate_consensus:
    input:
        lemur=rules.generate_consensus.output.fasta,
        others=consensus_dir / "non_lemur.consensus.nanopore.fa",
        masked_ref=rules.mask_h37rv.output.masked_genome,
    output:
        fasta=consensus_dir / "consensus.nanopore.fa",
    threads: 1
    resources:
        mem_mb=int(GB),
    log:
        rule_log_dir / "aggregate_consensus.log",
    shell:
        "awk 1 {input.lemur} {input.others} {input.masked_ref} > {output.fasta} 2> {log}"


rule lemur_snp_distance:
    input:
        fasta=rules.aggregate_consensus.output.fasta,
    output:
        matrix=distance_dir / "matrix.nanopore.csv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(GB) * attempt,
    container:
        containers["snp-dists"]
    log:
        rule_log_dir / "nanopore_snp_distance.log",
    params:
        options="-c",
    shell:
        """
        snp-dists {params.options} {input.fasta} 2> {log} > {output.matrix}
        """


rule plot_lemur_distance_matrix:
    input:
        matrix=rules.lemur_snp_distance.output.matrix,
    output:
        plot=report(
            distance_dir / "heatmap.nanopore.html",
            category="Distance",
            caption=captions["distance_matrix"],
        ),
    params:
        script=scripts["plot_distance_matrix"],
        options=" ".join(["--delim ,", "--title 'Pairwise distance'",]),
    threads: 1
    resources:
        mem_mb=int(GB),
    container:
        containers["conda"]
    conda:
        envs["plot_distance_matrix"]
    shell:
        """
        python {params.script} {params.options} -i {input.matrix} -o {output.plot}
        """


rule phylogenetic_tree:
    input:
        alignment=rules.aggregate_consensus.output.fasta,
    output:
        tree=phylo_dir / "lemur.tree",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(16 * GB) * attempt,
    container:
        containers["fasttree"]
    log:
        rule_log_dir / "phylogenetic_tree.log",
    params:
        options="-gtr -nt -fastest",
    shell:
        "FastTree {params.options} {input.alignment} > {output.tree} 2> {log}"
