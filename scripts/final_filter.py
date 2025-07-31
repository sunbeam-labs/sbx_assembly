from itertools import groupby
from typing import Iterator, TextIO, Tuple


def parse_fasta(f: TextIO) -> Iterator[Tuple[str, str]]:
    faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq_str = "".join(s.strip() for s in faiter.__next__())

        yield (header_str, seq_str)


def write_fasta(record: Tuple[str, str], f: TextIO) -> None:
    f.write(f">{record[0]}\n")
    f.write(f"{record[1]}\n")


with open(snakemake.log[0], "w") as log:
    with open(snakemake.input[0]) as f_in, open(snakemake.output[0], "w") as f_out:
        for record in parse_fasta(f_in):
            if len(record[1]) >= snakemake.params.len:
                write_fasta(record, f_out)
            else:
                log.write(f"Filtered {record[0]} for being too short")
