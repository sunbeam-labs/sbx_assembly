import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile


@pytest.fixture
def setup():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")
    hosts_fp = os.path.abspath(".tests/data/hosts/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    output_fp = os.path.join(project_dir, "sunbeam_output")

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "--target_list",
                "all_assembly",
                "all_annotate",
                "--directory",
                temp_dir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
        sys.exit(e)

    try:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
    except FileExistsError:
        pass

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield output_fp, benchmarks_fp


def test_full_run_assembly(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    final_contigs_fp = os.path.join(output_fp, "assembly/contigs/TEST-contigs.fa")
    genes_fp = os.path.join(output_fp, "annotation/genes/prodigal")

    # Check output
    assert os.path.exists(final_contigs_fp)
    assert os.stat(final_contigs_fp).st_size > 0
    for ext in ["_nucl.fa", "_prot.fa", ".gff"]:
        assert os.path.exists(os.path.join(genes_fp, f"TEST_genes{ext}"))


def test_full_run_annotation(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    # Check output
    all_samples_fp = os.path.join(output_fp, "annotation/all_samples.tsv")
    assert os.path.exists(all_samples_fp)
    assert os.stat(all_samples_fp).st_size > 0
