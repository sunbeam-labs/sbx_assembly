from sunbeam.bfx.parse import parse_fasta, write_fasta

with open(snakemake.log[0], "w") as log:
    with open(snakemake.input[0]) as f_in, open(snakemake.output[0], "w") as f_out:
        for record in parse_fasta(f_in):
            if len(record[1]) >= snakemake.params.len:
                write_fasta(record, f_out)
            else:
                log.write(f"Filtered {record[0]} for being too short")
