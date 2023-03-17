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
        "sbx_coverage.yml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.contig} {input.reads} > {output} 2>&1 | tee {log}
        """


rule contigs_sort:
    input:
        ASSEMBLY_FP / "contigs" / "minimap2" / "{sample}.sam",
    output:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.sam",
    benchmark:
        BENCHMARK_FP / "contigs_sort_{sample}.tsv"
    log:
        LOG_FP / "contigs_sort_{sample}.log",
    threads: Cfg["sbx_coverage"]["threads"]
    # conda:
    #    "sbx_coverage.yml"
    wrapper:
        "v1.23.5/bio/samtools/sort"


# shell:
#    """
#    samtools sort -@ {threads} -o {output.bam} {input} 2>&1 | tee {log}
#    samtools index {output.bam} {output.bai} 2>&1 | tee {log}
#    """


rule contigs_index:
    input:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.sam",
    output:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.sam.bai",
    benchmark:
        BENCHMARK_FP / "contigs_index_{sample}.tsv"
    log:
        LOG_FP / "contigs_index_{sample}.log",
    wrapper:
        "v1.23.5/bio/samtools/index"


rule mapping_depth:
    input:
        ASSEMBLY_FP / "contigs" / "samtools" / "{sample}.sorted.sam",
    output:
        ASSEMBLY_FP / "contigs" / "coverage" / "{sample}.depth",
    benchmark:
        BENCHMARK_FP / "mapping_depth_{sample}.tsv"
    log:
        LOG_FP / "mapping_depth_{sample}.log",
    # conda:
    #    "sbx_coverage.yml"
    wrapper:
        "v1.23.5/bio/samtools/depth"


# shell:
#    """
#    samtools depth -aa {input.bam} -o {output.depth} 2>&1 | tee {log}
#    """


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
        "sbx_coverage.yml"
    script:
        "scripts/assembly/get_coverage.py"


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
