import os
import subprocess
import pytest

# import numpy as np


@pytest.fixture()
def run_wrstart(get_chain, submit_options):
    if submit_options is None:
        submit_options = ""
    # we actually need to wait for submit to finish
    # wrstart finishes first then wruniq
    # for wruniq, in wruniq1.log there appears a line
    # ""
    submit_path = "${POWR_WORK}/proc.dir/submit.com wrstart1"
    submit_command = submit_path + submit_options
    # run wrstart

    try:
        subprocess.run(
            submit_command,
            shell=True,
            check=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        print(error.stderr)
        print(error.stdout)
        assert False, "CalledProcessError error"

    os.system("ls ${POWR_WORK}")
    os.system("ls ${POWR_WORK}/wrdata1")
    os.system("ls ${POWR_WORK}/output")
    os.system("cat ${POWR_WORK}/output/*.cpr")
    yield "ran powr full cycle"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    return "Cleaned powr tmp data"


@pytest.mark.parametrize("submit_options", ["", " nonopt"])
def test_run_wrstart(run_wrstart):
    print("done")
