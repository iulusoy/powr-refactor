import os
import subprocess
from pathlib import Path
import pytest
# import numpy as np


# get the path of this file and export variables accordingly
@pytest.fixture(scope="module")
def set_vars():
    filedir = Path(__file__).parents[0]
    powrdir = filedir.parents[0] / "powr"
    os.environ["POWR_WORK"] = powrdir.as_posix()
    os.environ["POWREXEPATH"] = (powrdir / "exe.dir").as_posix()
    print("powr work in set_vars")
    os.system("echo ${POWR_WORK}")
    return powrdir


@pytest.fixture(scope="module")
def set_aliases(set_vars):
    setup_file = "bash_setup"  # this could be different on MacOS
    # the aliases setup currently does not work as the aliases
    # are destryed between subprocess sessions
    # either we set them in the python script or
    # run explicitly with passing the env as dict
    full_file_path = "${POWR_WORK}/proc.dir/" + setup_file
    subprocess.run(
        full_file_path,
        shell=True,
        check=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )


@pytest.fixture(scope="module")
def get_chain(set_aliases):
    makechain_command = "${POWR_WORK}/proc.dir/makechain.com 1"
    subprocess.run(
        makechain_command,
        shell=True,
        check=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )
    os.system("ls ${POWR_WORK}")
    yield "Created chain 1"
    # teardown directories
    # we need access to ${POWR_WORK} so shutil will not work
    # why does sourcing powrconfig in the script not work?
    os.system("rm -rf ${POWR_WORK}/scratch")
    os.system("rm -rf ${POWR_WORK}/output")
    os.system("rm -rf ${POWR_WORK}/wrdata1")
    os.system("rm -rf ${POWR_WORK}/wrjobs")
    return "Cleaned working dirs"


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