with open(snakemake.output[0], "w") as out:
    out.writelines(open(snakemake.input[0]).readlines())
    for infile in snakemake.input[1:]:
        out.writelines(l for i, l in enumerate(open(infile)) if i > 0)
