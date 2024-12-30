import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_random_forest(pytestconfig):

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/random_forest/data")
        expected_path = PurePosixPath("tests/random_forest/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "workflow-results/EXPLANA-install-test/SELECTED-FEATURES-original/original-boruta-important.txt",
            file=sys.stderr
        )

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "workflow-results/EXPLANA-install-test/SELECTED-FEATURES-original/original-boruta-important.txt",
            "-f", 
            "-j1",
            "--target-files-omit-workdir-adjustment",
            "--configfile",
            f"{pytestconfig.rootdir}/config/config-install-test.yaml",
            "--directory",
            workdir,
        ])

        df = pd.read_csv(
            workdir / "workflow-results/EXPLANA-install-test/SELECTED-FEATURES-original/original-boruta-important.txt",
            sep="\t"
        )
        df_exp = pd.read_csv(
            expected_path / "workflow-results/EXPLANA-install-test/SELECTED-FEATURES-original/original-boruta-important.txt",
            sep="\t"
        )
        assert set(df.important_features.iloc[0:6]) == set(df_exp.important_features.iloc[0:6])
        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        #common.OutputChecker(data_path, expected_path, workdir).check()
