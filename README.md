<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_assembly

<!-- Begin badges -->
[![Tests](https://github.com/sunbeam-labs/sbx_assembly/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_assembly/actions/workflows/tests.yml)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_assembly)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_assembly/)
<!-- End badges -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for assembly of contigs using Megahit, gene annotation using Prodigal, and annotation using [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [Diamond](https://github.com/bbuchfink/diamond). It can also map reads to contigs and calculat per-base coverage using [Minimap2](https://github.com/lh3/minimap2) and [samtools](https://github.com/samtools/samtools).

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_assembly.git

## Usage

To generate contigs, annotations and coverage, create a project and use the `all_assembly`, `all_annotation`, and `all_coverage` targets:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    sunbeam run --profile /path/to/project/ --target_list all_assembly all_annotate

N.B. For sunbeam versions <4 the last command will be something like `sunbeam run --configfile /path/to/project/sunbeam_config.yml all_assembly` and `sunbeam run --configfile /path/to/project/sunbeam_config.yml all_annotation`.

## Configuration

Assembly
  - min_length: Is the minimum contig length to accept in the final filtering
Annotation
  - threads: Is the number of threads to use running blast/diamond
  - circular_kmin: Is the minimum length of overlap to check for circular contigs with
  - circular_kmax: Is the maximum length of overlap to check for circular contigs with
  - circular_min_len: Is the minimum sequence length to consider checking for circularity
Coverage
  - threads: Is the number of threads to use while running `minimap2` and `samtools sort`

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_assembly.git extensions/sbx_assembly

and then include it in the config for any given project with:

    cat extensions/sbx_assembly/config.yml >> /path/to/project/sunbeam_config.yml
