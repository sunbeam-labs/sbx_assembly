# -*- mode: Snakemake -*-
#
# Contig building and other assembly rules
#
# Requires Megahit.


TARGET_ASSEMBLY = [
    expand(ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa", sample=Samples.keys()),
    expand(
        ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_{suffix}.fa",
        sample=Samples.keys(),
        suffix=["prot", "nucl"],
    ),
]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


rule all_assembly:
    """Build contigs for all samples."""
    input:
        TARGET_ASSEMBLY,


ruleorder: megahit_paired > megahit_unpaired


rule megahit_paired:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    benchmark:
        BENCHMARK_FP / "megahit_paired_{sample}.tsv"
    log:
        LOG_FP / "megahit_paired_{sample}.log",
    params:
        out_fp=str(ASSEMBLY_FP / "megahit" / "{sample}_asm"),
    threads: 4
    conda:
        "sbx_assembly.yml"
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        ## sometimes the error is due to lack of memory
        exitcode=0
        if [ -d {params.out_fp} ]
        then
            echo "Clearing previous megahit directory..." > {log}
            rm -rf {params.out_fp}
        fi
        megahit -t {threads} -1 {input.r1} -2 {input.r2} -o {params.out_fp} --continue 2>&1 {log} || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            touch {output}
            echo "Empty contigs" 2>&1 | tee {log}
        elif [ $exitcode -gt 1 ]
        then
            echo "Check your memory" 2>&1 | tee {log}
        fi
        """


rule megahit_unpaired:
    input:
        QC_FP / "decontam" / "{sample}_1.fastq.gz",
    output:
        ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    benchmark:
        BENCHMARK_FP / "megahit_unpaired_{sample}.tsv"
    log:
        LOG_FP / "megahit_unpaired_{sample}.log",
    params:
        out_fp=str(ASSEMBLY_FP / "megahit" / "{sample}_asm"),
    threads: 4
    conda:
        "sbx_assembly.yml"
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        ## sometimes the error is due to lack of memory
        exitcode=0
        megahit -t {threads} -r {input} -o {params.out_fp} -f --continue 2>&1 {log} || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            echo "Empty contigs"
            touch {output}
        elif [ $exitcode -gt 1 ]
        then
            echo "Check your memory"
        fi
        """


rule final_filter:
    input:
        ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    log:
        LOG_FP / "final_filter_{sample}.log",
    benchmark:
        BENCHMARK_FP / "final_filter_{sample}.tsv"
    params:
        len=Cfg["sbx_assembly"]["min_length"],
    script:
        "scripts/final_filter.py"


rule clean_assembly:
    input:
        M=ASSEMBLY_FP / "megahit",
    shell:
        """
        rm -rf {input.M} && echo "Cleanup assembly finished."
        """


rule prodigal:
    """Use Progial for coding genes predictions in contigs."""
    input:
        ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    output:
        gff=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes.gff",
        faa=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_prot.fa",
        fna=ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_nucl.fa",
    benchmark:
        BENCHMARK_FP / "prodigal_{sample}.tsv"
    log:
        LOG_FP / "prodigal_{sample}.log",
    conda:
        "sbx_assembly.yml"
    shell:
        """
        if [[ -s {input} ]]; then
          prodigal -i {input} -o {output.gff} \
          -a {output.faa} -d {output.fna} -p meta 2>&1 | tee {log}
        else
          touch {output.faa}
          touch {output.gff}
          touch {output.fna}
        fi
        """


rule _test_prodigal:
    input:
        expand(
            ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_{suffix}.fa",
            sample=Samples.keys(),
            suffix=["prot", "nucl"],
        ),