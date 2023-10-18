import os
import subprocess
from pathlib import Path
import pytest

# get the path of this file and export variables accordingly
@pytest.fixture
def set_vars():
    filedir = Path(__file__).parents[0]
    powrdir = filedir.parents[0] / "powr"
    os.environ["POWR_WORK"] = powrdir.as_posix()
    os.environ["POWREXEPATH"] = (powrdir / "exe.dir").as_posix()
    return powrdir

@pytest.fixture
def set_aliases(set_vars):
    setup_file = "bash_setup" # this could be different on MacOS
    full_file_path = "${POWR_WORK}/proc.dir/" + setup_file
    subprocess.run(full_file_path, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)

@pytest.fixture
def get_chain(set_aliases):
    makechain_command = "${POWR_WORK}/proc.dir/makechain.com 1"
    subprocess.run(makechain_command, shell=True, check=True, executable='/bin/bash', capture_output=True, text=True)
    yield "Created chain 1"
    # teardown directories
    os.system("rm -rf ${POWR_WORK}/scratch")
    os.system("rm -rf ${POWR_WORK}/output")
    os.system("rm -rf ${POWR_WORK}/wrdata1")
    os.system("rm -rf ${POWR_WORK}/wrjobs")
    return "Cleaned working dirs"


# test the correct set-up of folder structure for jobs
def test_makechain(set_vars, get_chain):
    assert set_vars.exists()

# compare the output from the colitest job
def test_read_model():
    pass