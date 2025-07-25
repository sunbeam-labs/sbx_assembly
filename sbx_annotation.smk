try:
    SBX_ASSEMBLY_VERSION = get_ext_version("sbx_assembly")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_ASSEMBLY_VERSION = "0.0.0"

try:
    SBX_KRAKEN_VERSION = get_ext_version("sbx_kraken")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_KRAKEN_VERSION = "0.0.0"

rule all_annotate:
    input:
        expand(ASSEMBLY_FP / "bakta" / "{sample}" / "{sample}.txt", sample=Samples),
        expand(ASSEMBLY_FP / "kraken" / "{sample}-raw.tsv", sample=Samples),

rule mag_bakta:
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    output:
        bakta=ASSEMBLY_FP / "bakta" / "{sample}" / "{sample}.txt",
    params:
        ref=Cfg["sbx_assembly"]["bakta_ref"],
    log:
        LOG_FP / "assembly_bakta_{sample}.log",
    benchmark:
        BENCHMARK_FP / "assembly_bakta_{sample}.tsv"
    conda:
        "envs/sbx_annotation.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            bakta --force --db {params.ref} \\
            --output $(dirname {output.bakta}) \\
            --prefix {wildcards.sample} \\
            --skip-plot {input.contigs} \\
            &> {log}
        else
            touch {output.bakta}
        fi
        """


rule kraken2_classify_mags:
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    output:
        raw=ASSEMBLY_FP / "kraken" / "{sample}-raw.tsv",
        report=ASSEMBLY_FP / "kraken" / "report" / "{sample}-taxa.tsv",
    benchmark:
        BENCHMARK_FP / "kraken2_assembly_{sample}.tsv"
    log:
        LOG_FP / "kraken2_assembly_{sample}.log",
    params:
        db=Cfg["sbx_assembly"]["kraken_db_fp"],
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_kraken:{SBX_KRAKEN_VERSION}"
    threads: 8
    shell:
        """
        if [ -s {input.contigs} ]; then
            kraken2 \
            --db {params.db} \
            --report {output.report} \
            --output {output.raw} \
            {input} \
            > {log} 2>&1
        else
            echo "Empty input, filling with dummy data" > {log}
            echo "100.00\t0\t0\tR\t1\troot" > {output.report}
            echo "C\tA\t1\t136|136\t1:102 |:| 1:102" > {output.raw}
        fi
        """
