import csv
from collections import Counter
from sunbeam.bfx.parse import parse_fasta
from . import blast_summary


def circular(seq: str, kmin: int, kmax: int, min_len: int) -> bool:
    """Determine if a sequence is circular.

    Checks for repeated k-mer at beginning and end of a sequence for a given
    range of values for k.
    """
    if len(seq) < min_len:
        return False
    # Short-circuit checking: returns True for the first kmer that matches
    return any([k for k in range(kmin, kmax + 1) if seq[0:k] == seq[len(seq) - k :]])


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
