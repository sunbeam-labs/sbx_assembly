try:
    SBX_ASSEMBLY_VERSION = get_ext_version("sbx_assembly")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_ASSEMBLY_VERSION = "0.0.0"


rule all_assembly:
    """Build contigs for all samples."""
    input:
        expand(ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa", sample=Samples.keys()),
        expand(
            ANNOTATION_FP / "genes" / "prodigal" / "{sample}_genes_{suffix}.fa",
            sample=Samples.keys(),
            suffix=["prot", "nucl"],
        ),


ruleorder: megahit_paired > megahit_unpaired


rule megahit_paired:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        final=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    benchmark:
        BENCHMARK_FP / "megahit_paired_{sample}.tsv"
    log:
        LOG_FP / "megahit_paired_{sample}.log",
    threads: 4
    conda:
        "envs/sbx_assembly.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-assembly"
    shell:
        """
        # Turn off bash strict mode
        set +o pipefail

        OUT_DIR=$(dirname {output.final})
        if [ -d $OUT_DIR ]
        then
            echo "Clearing previous megahit directory..." > {log}
            rm -rf $OUT_DIR
        fi

        exitcode=0
        megahit -t {threads} -1 {input.r1} -2 {input.r2} -o $OUT_DIR --continue >> {log} 2>&1 || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            touch {output}
            echo "Empty contigs" >> {log}
        elif [ $exitcode -gt 1 ]
        then
            echo "Something went wrong (maybe check the memory usage)" >> {log}
        fi
        """


rule megahit_unpaired:
    input:
        QC_FP / "decontam" / "{sample}_1.fastq.gz",
    output:
        final=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    benchmark:
        BENCHMARK_FP / "megahit_unpaired_{sample}.tsv"
    log:
        LOG_FP / "megahit_unpaired_{sample}.log",
    threads: 4
    conda:
        "envs/sbx_assembly.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-assembly"
    shell:
        """
        # Turn off bash strict mode
        set +o pipefail

        OUT_DIR=$(dirname {output.final})
        if [ -d $OUT_DIR ]
        then
            echo "Clearing previous megahit directory..." > {log}
            rm -rf $OUT_DIR
        fi

        exitcode=0
        megahit -t {threads} -r {input} -o $OUT_DIR -f --continue >> {log} 2>&1 || exitcode=$?

        if [ $exitcode -eq 255 ]
        then
            echo "Empty contigs" >> {log}
            touch {output}
        elif [ $exitcode -gt 1 ]
        then
            echo "Something went wrong (maybe check the memory usage)" >> {log}
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
        "envs/sbx_assembly.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-assembly"
    shell:
        """
        if [[ -s {input} ]]; then
          prodigal -i {input} -o {output.gff} \
          -a {output.faa} -d {output.fna} -p meta > {log} 2>&1
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
