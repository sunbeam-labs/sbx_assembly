import os
import pytest
import shutil
import subprocess as sp
import tempfile
from pathlib import Path


@pytest.fixture
def setup():
    temp_dir = Path(tempfile.mkdtemp())

    reads_fp = Path(".tests/data/reads/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()
    project_dir = temp_dir / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"qc: {{host_fp: {hosts_fp}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_annotate",
            "all_assembly",
            "all_coverage",
            "--directory",
            temp_dir,
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    if os.getenv("GITHUB_ACTIONS") == "true":
        try:
            shutil.copytree(log_fp, "logs/")
            shutil.copytree(stats_fp, "stats/")
        except FileNotFoundError:
            print("No logs or stats directory found.")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam run failed with error: {proc.stderr}"

    ### Assembly
    lfinal_contigs_fp = output_fp / "assembly" / "contigs" / "LONG-contigs.fa"
    sfinal_contigs_fp = output_fp / "assembly" / "contigs" / "SHORT-contigs.fa"
    genes_fp = output_fp / "annotation" / "genes" / "prodigal"

    assert lfinal_contigs_fp.exists()
    assert lfinal_contigs_fp.stat().st_size > 0
    assert sfinal_contigs_fp.exists()
    for ext in ["_nucl.fa", "_prot.fa", ".gff"]:
        assert (genes_fp / f"LONG_genes{ext}").exists()
        assert (genes_fp / f"SHORT_genes{ext}").exists()

    ### Annotation
    all_samples_fp = output_fp / "annotation" / "all_samples.tsv"
    blastn_fp = output_fp / "annotation" / "blastn" / "bacteria" / "contig" / "LONG.btf"
    blastp_fp = output_fp / "annotation" / "blastp" / "prot" / "prodigal" / "LONG.btf"
    blastx_fp = output_fp / "annotation" / "blastx" / "prot" / "prodigal" / "LONG.btf"

    assert all_samples_fp.exists()
    assert all_samples_fp.stat().st_size > 0

    ### Coverage
    contigs_coverage_fp = output_fp / "assembly" / "contigs_coverage.txt"

    assert contigs_coverage_fp.exists()
    assert contigs_coverage_fp.stat().st_size > 0
