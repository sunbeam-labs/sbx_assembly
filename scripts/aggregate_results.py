import csv
from collections import Counter
from pathlib import Path
from sunbeamlib import circular
from sunbeamlib.parse import parse_fasta
from xml.etree.ElementTree import ParseError


BLAST6_DEFAULTS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart," "qend," "sstart," "send",
    "evalue",
    "bitscore",
]


def parse_blast6(f, outfmt=BLAST6_DEFAULTS):
    for line in f.readlines():
        vals = line.strip().split("\t")
        if len(outfmt) == len(vals):
            yield dict(zip(outfmt, vals))


def blast_summary(blast_files):
    """Summarize BLAST results from an set of BLAST output files."""
    for infile in blast_files:
        sample = Path(infile).stem
        try:
            with open(infile) as f:
                for result in parse_blast6(f):
                    if len(result.hits) > 0:
                        yield {
                            "sample": sample,
                            "query": result.id,
                            "hit": result.hits[0].id,
                        }
        except ParseError:
            print("Skipping empty/malformed %s" % infile)
            continue


def blast_contig_summary(blast_files):
    return {r["query"]: r["hit"] for r in blast_summary(blast_files)}


def blast_hits(blast_files):
    return Counter(c["query"] for c in blast_summary(blast_files))


with open(snakemake.log[0], "w") as log:
    log.write("Starting aggregation...\n")
    with open(snakemake.input.contigs) as f:
        contigs = {r[0]: r[1] for r in parse_fasta(f)}
    # Separate each set of result files by the database it was blasted against
    contig_results = {
        db: blast_contig_summary(f for f in snakemake.input.contig_results if db in f)
        for db in snakemake.params.nucl
    }
    # We only care about the number of hits for the protein database
    gene_hits = {
        db: blast_hits(f for f in snakemake.input.gene_results if db in f)
        for db in snakemake.params.prot
    }

    with open(snakemake.output[0], "w") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["sample", "contig", "length", "circular"]
            + snakemake.params.dbs,
            delimiter="\t",
        )
        writer.writeheader()
        for contig, contig_seq in contigs.items():
            is_circular = circular(
                contig_seq,
                snakemake.config["sbx_annotation"]["circular_kmin"],
                snakemake.config["sbx_annotation"]["circular_kmax"],
                snakemake.config["sbx_annotation"]["circular_min_len"],
            )
            results = {
                "sample": snakemake.wildcards.sample,
                "contig": contig,
                "length": len(contig_seq),
                "circular": is_circular,
            }
            for db in snakemake.params.nucl:
                results[db] = contig_results[db].get(contig, "NA")
            # Report the number of hits of each contig/gene for each prot. db
            # Genes are reported from prodigal as contig_1,contig_2, etc so
            # aggregate all hits for all of a contig's genes together using sum
            for db in snakemake.params.prot:
                results[db] = sum(
                    gene_hits[db][gene] for gene in gene_hits[db] if contig in gene
                )
            writer.writerow(results)
