import csv
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))
from scripts import blast_summary


with open(snakemake.output[0], "w") as out:
    writer = csv.DictWriter(out, fieldnames=["sample", "query", "hit"], delimiter="\t")
    writer.writeheader()
    list(writer.writerow(result) for result in blast_summary(snakemake.input))
