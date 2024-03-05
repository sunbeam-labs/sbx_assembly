FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_assembly_env

COPY envs/sbx_annotation.yml ./

# Install environment
RUN mamba env create --file sbx_annotation.yml --name sbx_annotation

ENV PATH="/opt/conda/envs/sbx_annotation/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_annotation", "/bin/bash", "-c"]

# Run
CMD "bash"