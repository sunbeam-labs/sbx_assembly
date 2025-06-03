try:
    SBX_ASSEMBLY_VERSION = get_ext_version("sbx_assembly")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_ASSEMBLY_VERSION = "0.0.0"


Blastdbs = {"nucl": {}, "prot": {}}
try:
    Blastdbs["nucl"] = {
        db: str(Path(Cfg["blastdbs"]["root_fp"]) / path)
        for db, path in Cfg["blastdbs"].get("nucleotide", {}).items()
    }
    Blastdbs["prot"] = {
        db: str(Path(Cfg["blastdbs"]["root_fp"]) / path)
        for db, path in Cfg["blastdbs"].get("protein", {}).items()
    }
except KeyError:
    pass


rule all_annotate:
    input:
        ANNOTATION_FP / "all_samples.tsv",


rule build_diamond_db:
    """Use diamond makedb to create any necessary db indeces that don't exist."""
    input:
        [Blastdbs["prot"][db] for db in Blastdbs["prot"]],
    output:
        [Blastdbs["prot"][db] + ".dmnd" for db in Blastdbs["prot"]],
    benchmark:
        BENCHMARK_FP / "build_diamond_db.tsv"
    log:
        LOG_FP / "build_diamond_db.log",
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    shell:
        """
        diamond makedb --in {input} -d {input} > {log} 2>&1
        """


rule run_blastn:
    """Run BLASTn against a given database and write the results to blast tabular format."""
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
    output:
        ANNOTATION_FP / "blastn" / "{db}" / "{contigs}" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_blastn_{db}_{contigs}_{sample}.tsv"
    log:
        LOG_FP / "run_blastn_{db}_{contigs}_{sample}.log",
    params:
        db=lambda wildcard: Blastdbs["nucl"][wildcard.db],
    threads: Cfg["sbx_annotation"]["threads"]
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    shell:
        """
        blastn \
        -query {input.contigs} \
        -db {params.db} \
        -outfmt 7 \
        -num_threads {threads} \
        -evalue 1e-10 \
        -max_target_seqs 5000 \
        -out {output} \
        2>&1 | tee {log}
        """


rule run_diamond_blastp:
    """Run diamond blastp on translated genes against a target db and write to blast tabular format."""
    input:
        genes=ANNOTATION_FP / "genes" / "{orf_finder}" / "{sample}_genes_prot.fa",
        indexes=rules.build_diamond_db.output,
    output:
        ANNOTATION_FP / "blastp" / "{db}" / "{orf_finder}" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_diamond_blastp_{db}_{orf_finder}_{sample}.tsv"
    log:
        LOG_FP / "run_diamond_blastp_{db}_{orf_finder}_{sample}.log",
    threads: Cfg["sbx_annotation"]["threads"]
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    shell:
        """
        if [ -s {input.genes} ]; then
            diamond blastp \
            -q {input.genes} \
            --db {input.indexes} \
            --outfmt 6 \
            --threads {threads} \
            --evalue 1e-10 \
            --max-target-seqs 2475 \
            --out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


rule run_diamond_blastx:
    """Run diamond blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        genes=ANNOTATION_FP / "genes" / "{orf_finder}" / "{sample}_genes_nucl.fa",
        indexes=rules.build_diamond_db.output,
    output:
        ANNOTATION_FP / "blastx" / "{db}" / "{orf_finder}" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_diamond_blastx_{db}_{orf_finder}_{sample}.tsv"
    log:
        LOG_FP / "run_diamond_blastx_{db}_{orf_finder}_{sample}.log",
    threads: Cfg["sbx_annotation"]["threads"]
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    shell:
        """
        if [ -s {input.genes} ]; then
            diamond blastx \
            -q {input.genes} \
            --db {input.indexes} \
            --outfmt 6 \
            --threads {threads} \
            --evalue 1e-10 \
            --max-target-seqs 2475 \
            --out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


rule blast_report:
    """Create a summary of results from a BLAST call."""
    input:
        expand(
            ANNOTATION_FP / "{{blast_prog}}" / "{{db}}" / "{{query}}" / "{sample}.btf",
            sample=Samples.keys(),
        ),
    output:
        ANNOTATION_FP / "{blast_prog}" / "{db}" / "{query}" / "report.tsv",
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    script:
        "scripts/blast_report.py"


rule _test_blastpx:
    input:
        expand(
            ANNOTATION_FP / "{blastpx}" / "card" / "prodigal" / "{sample}.btf",
            blastpx=["blastx", "blastp"],
            sample=Samples.keys(),
        ),


rule _test_blastpx_report:
    input:
        expand(
            ANNOTATION_FP / "{blastpx}" / "card" / "prodigal" / "report.tsv",
            blastpx=["blastx", "blastp"],
        ),


rule aggregate_results:
    input:
        contigs=ASSEMBLY_FP / "contigs" / "{sample}-contigs.fa",
        contig_results=expand(
            ANNOTATION_FP / "blastn" / "{db}" / "contig" / "{{sample}}.btf",
            db=Blastdbs["nucl"],
        ),
        gene_results=expand(
            ANNOTATION_FP / "{blastpx}" / "{db}" / "{orf_finder}" / "{{sample}}.btf",
            blastpx=["blastp", "blastx"],
            db=Blastdbs["prot"],
            orf_finder=["prodigal"],
        ),
    output:
        ANNOTATION_FP / "summary" / "{sample}.tsv",
    benchmark:
        BENCHMARK_FP / "aggregate_results_{sample}.tsv"
    log:
        LOG_FP / "aggregate_results_{sample}.log",
    params:
        dbs=list(Blastdbs["nucl"].keys()) + list(Blastdbs["prot"].keys()),
        nucl=Blastdbs["nucl"],
        prot=Blastdbs["prot"],
    conda:
        "envs/sbx_annotation.yml"
    container:
        f"docker://sunbeamlabs/sbx_assembly:{SBX_ASSEMBLY_VERSION}-annotation"
    script:
        "scripts/aggregate_results.py"


rule aggregate_all:
    input:
        expand(ANNOTATION_FP / "summary" / "{sample}.tsv", sample=Samples.keys()),
    output:
        ANNOTATION_FP / "all_samples.tsv",
    script:
        "scripts/aggregate_all.py"
