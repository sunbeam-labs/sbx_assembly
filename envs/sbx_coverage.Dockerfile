FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_assembly_env

COPY envs/sbx_coverage.yml ./

# Install environment
RUN mamba env create --file sbx_coverage.yml --name sbx_coverage

ENV PATH="/opt/conda/envs/sbx_coverage/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_coverage", "/bin/bash", "-c"]

# Run
CMD "bash"