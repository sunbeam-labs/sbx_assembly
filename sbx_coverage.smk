# -*- mode: Snakemake -*-
#
# Map reads to contigs and calculate per base coverage
#
# Requires Minimap2 and samtools.

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


rule all_coverage:
    input:
        ASSEMBLY_FP / "contigs_coverage.txt",


rule _contigs_mapping:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
            sample=Samples.keys(),
        ),


rule _all_coverage:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
            sample=Samples.keys(),
        ),


rule minimap_alignment:
    input:
        contig=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        temp(ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sam"),
    benchmark:
        BENCHMARK_FP / "minimap_alignment_{sample}.tsv"
    log:
        LOG_FP / "minimap_alignment_{sample}.log",
    threads: Cfg["sbx_coverage"]["threads"]
    conda:
        "envs/sbx_coverage.yml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.contig} {input.reads} 1> {output} 2> {log}
        """


rule contigs_sort:
    input:
        ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sam",
    output:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.bam",
    benchmark:
        BENCHMARK_FP / "contigs_sort_{sample}.tsv"
    log:
        LOG_FP / "contigs_sort_{sample}.log",
    threads: Cfg["sbx_coverage"]["threads"]
    conda:
        "envs/sbx_coverage.yml"
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2>&1 | tee {log}
        """


rule mapping_depth:
    input:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.bam",
    output:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
    benchmark:
        BENCHMARK_FP / "mapping_depth_{sample}.tsv"
    log:
        LOG_FP / "mapping_depth_{sample}.log",
    conda:
        "envs/sbx_coverage.yml"
    shell:
        """
        samtools depth -aa {input} 1> {output} 2> {log}
        """


rule get_coverage:
    input:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
    output:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
    benchmark:
        BENCHMARK_FP / "get_coverage_{sample}.tsv"
    log:
        LOG_FP / "get_coverage_{sample}.log",
    conda:
        "envs/sbx_coverage.yml"
    script:
        "scripts/get_coverage.py"


rule summarize_coverage:
    input:
        expand(
            ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.csv",
            sample=Samples.keys(),
        ),
    output:
        ASSEMBLY_FP / "contigs_coverage.txt",
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"
