import os
import subprocess
import pytest

# import numpy as np


@pytest.fixture(scope="module")
def run_wrstart(get_chain):
    # run wrstart
    wrstart_command = "${POWR_WORK}/proc.dir/submit.com wrstart1"
    temp = subprocess.run(
        wrstart_command,
        shell=True,
        check=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )
    print(temp.stdout)
    print(temp.stderr)
    os.system("ls ${POWR_WORK}")
    os.system("ls ${POWR_WORK}/wrdata1")
    os.system("ls ${POWR_WORK}/output")
    os.system("cat ${POWR_WORK}/output/*.cpr")

    yield "ran powr full cycle"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    return "Cleaned powr tmp data"


def test_run_wrstart(run_wrstart):
    print("done")
