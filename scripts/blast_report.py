import csv
from pathlib import Path
from sunbeamlib.parse import parse_blast6
from xml.etree.ElementTree import ParseError


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


with open(snakemake.output[0], "w") as out:
    writer = csv.DictWriter(out, fieldnames=["sample", "query", "hit"], delimiter="\t")
    writer.writeheader()
    list(writer.writerow(result) for result in blast_summary(snakemake.input))
