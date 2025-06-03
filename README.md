<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_assembly

<!-- Begin badges -->
[![Tests](https://github.com/sunbeam-labs/sbx_assembly/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_assembly/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_assembly)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_assembly/)
<!-- End badges -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for assembly of contigs using Megahit, gene annotation using Prodigal, and annotation using [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [Diamond](https://github.com/bbuchfink/diamond). It can also map reads to contigs and calculat per-base coverage using [Minimap2](https://github.com/lh3/minimap2) and [samtools](https://github.com/samtools/samtools).

## Configuration

sbx_assembly
  - min_length: Is the minimum contig length to accept in the final filtering
abx_annotation
  - threads: Is the number of threads to use running blast/diamond
  - circular_kmin: Is the minimum length of overlap to check for circular contigs with
  - circular_kmax: Is the maximum length of overlap to check for circular contigs with
  - circular_min_len: Is the minimum sequence length to consider checking for circularity
sbx_coverage
  - threads: Is the number of threads to use while running `minimap2` and `samtools sort`

## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).