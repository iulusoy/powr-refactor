import pytest
import os
import subprocess
from pathlib import Path


# get the path of this file and export variables accordingly
@pytest.fixture(scope="session")
def set_vars():
    filedir = Path(__file__).parents[0]
    powrdir = filedir.parents[0] / "powr"
    os.environ["POWR_WORK"] = powrdir.as_posix()
    os.environ["POWREXEPATH"] = (powrdir / "exe_dev.dir").as_posix()
    # create the tmp_2day folder if not exists
    tmp_2day = powrdir / "tmp_2day"
    if not os.path.exists(tmp_2day):
        os.mkdir(tmp_2day)
    return powrdir


# inject path into powrconfig
@pytest.fixture(scope="session")
def inject_path(set_vars):
    powrconfig_file = set_vars / "powrconfig"
    with open(powrconfig_file, "r") as f:
        powrconfig = f.read()
    # now insert path in correct spot
    temp = powrconfig.split("export POWR_WORK=")
    powrconfig = temp[0] + "export POWR_WORK=" + set_vars.as_posix() + temp[1]
    powrconfig_file = set_vars / ".powrconfig"
    with open(powrconfig_file, "w") as f:
        f.write(powrconfig)
    os.system("cp ${POWR_WORK}/.powrconfig ${HOME}/.powrconfig")


@pytest.fixture(scope="session")
def get_chain(inject_path):

    makechain_command = "${POWR_WORK}/proc.dir/makechain.com 1 -force"

    try:
        subprocess.run(
            makechain_command,
            shell=True,
            check=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        print(error.stderr)
        print(error.stdout)

    yield "Created chain 1"
    # teardown directories
    # we need access to ${POWR_WORK} so shutil will not work
    # why does sourcing powrconfig in the script not work?
    os.system("rm -rf ${POWR_WORK}/scratch")
    os.system("rm -rf ${POWR_WORK}/output")
    os.system("rm -rf ${POWR_WORK}/wrdata1")
    os.system("rm -rf ${POWR_WORK}/wrjobs")
    return "Cleaned working dirs"
