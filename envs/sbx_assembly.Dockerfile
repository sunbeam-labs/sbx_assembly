FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_assembly_env

COPY envs/sbx_assembly.yml ./

# Install environment
RUN mamba env create --file sbx_assembly.yml --name sbx_assembly

ENV PATH="/opt/conda/envs/sbx_assembly/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_assembly", "/bin/bash", "-c"]

# Run
CMD "bash"