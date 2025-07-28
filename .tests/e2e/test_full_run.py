import os
import pytest
import shutil
import subprocess as sp
from pathlib import Path


@pytest.fixture
def setup(tmp_path):
    reads_fp = Path(".tests/data/reads/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()
    kraken_db_fp = Path(".tests/data/kraken_db/").resolve()
    project_dir = tmp_path / "project"

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

    config_str = f"sbx_assembly: {{kraken_db_fp: {kraken_db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmp_path, project_dir


@pytest.fixture
def run_sunbeam(setup):
    tmp_path, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            # We don't have a good way of testing Bakta yet
            "_test_kraken",
            "all_assembly",
            "all_coverage",
            "--directory",
            tmp_path,
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

    assert lfinal_contigs_fp.exists()
    assert lfinal_contigs_fp.stat().st_size > 0
    assert sfinal_contigs_fp.exists()

    ### Annotation
    lkraken_report_fp = output_fp / "assembly" / "kraken" / "report" / "LONG-taxa.tsv"

    assert lkraken_report_fp.exists()
    assert lkraken_report_fp.stat().st_size > 0

    ### Coverage
    contigs_coverage_fp = output_fp / "assembly" / "coverage" / "contigs_coverage.txt"

    assert contigs_coverage_fp.exists()
    assert contigs_coverage_fp.stat().st_size > 0
