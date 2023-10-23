import os
import subprocess
from pathlib import Path
import pytest


# get the path of this file and export variables accordingly
@pytest.fixture(scope="module")
def set_vars():
    filedir = Path(__file__).parents[0]
    powrdir = filedir.parents[0] / "powr"
    os.environ["POWR_WORK"] = powrdir.as_posix()
    os.environ["POWREXEPATH"] = (powrdir / "exe.dir").as_posix()
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
    yield "Created chain 1"
    # teardown directories
    # we need access to ${POWR_WORK} so shutil will not work
    os.system("rm -rf ${POWR_WORK}/scratch")
    os.system("rm -rf ${POWR_WORK}/output")
    os.system("rm -rf ${POWR_WORK}/wrdata1")
    os.system("rm -rf ${POWR_WORK}/wrjobs")
    return "Cleaned working dirs"


@pytest.fixture(scope="module")
def run_colitest(get_chain):
    # run colitest
    colitest_command = "${POWR_WORK}/wrjobs/colitest1"
    subprocess.run(
        colitest_command,
        shell=True,
        check=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )
    yield "ran colitest"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    return "Cleaned colitest tmp data"


# test the correct set-up of folder structure for jobs
def test_makechain(set_vars, get_chain):
    assert set_vars.exists()
    # now check that the scratch, output, wrdata1, wrjobs
    # dirs have been generated
    scratch_dir = set_vars / "scratch"
    output_dir = set_vars / "output"
    wrdata1_dir = set_vars / "wrdata1"
    wrjobs_dir = set_vars / "wrjobs"
    assert scratch_dir.exists()
    assert output_dir.exists()
    assert wrdata1_dir.exists()
    assert wrjobs_dir.exists()
    # now check that the correct input files have been copied into the dirs
    scratch_content = [
        "formal1",
        "modify1",
        "newdatom1",
        "newformal_cards1",
        "njn1",
        "steal1",
        "wrstart1",
        "wruniq1",
    ]
    output_content = []
    wrdata1_content = [
        "CARDS",
        "DATOM",
        "FGRID",
        "FORMAL_CARDS",
        "MODEL",
        "NEWDATOM_INPUT",
        "NEWFORMAL_CARDS_INPUT",
    ]
    wrjobs_content = [
        "coli_test",
        "formal_wrh_gen",
        "formal_wrh_xxl",
        "newformal_cards1",
        "set_repeat1",
        "tmphosts",
        "wrstart_wrh_hydro",
        "wruniq1",
        "wruniq_wrh_merged",
        "colitest1",
        "formal_wrh_hydro",
        "modify1",
        "newformal_cards_gen",
        "set_steal1",
        "wrstart1",
        "wrstart_wrh_merged",
        "wruniq_wrh_dev",
        "wruniq_wrh_vd20",
        "formal1",
        "formal_wrh_merged",
        "newdatom1",
        "njn1",
        "steal1",
        "wrstart_wrh_dev",
        "wrstart_wrh_vd20",
        "wruniq_wrh_gen",
        "wruniq_wrh_xxl",
        "formal_wrh_dev",
        "formal_wrh_vd20",
        "newdatom_gen",
        "njn_wrh_gen",
        "steal1_backup",
        "wrstart_wrh_gen",
        "wrstart_wrh_xxl",
        "wruniq_wrh_hydro",
    ]
    temp = [i.name for i in scratch_dir.iterdir()]
    assert sorted(temp) == sorted(scratch_content)
    temp = [i.name for i in output_dir.iterdir()]
    assert sorted(temp) == sorted(output_content)
    temp = [i.name for i in wrdata1_dir.iterdir()]
    assert sorted(temp) == sorted(wrdata1_content)
    temp = [i.name for i in wrjobs_dir.iterdir()]
    assert sorted(temp) == sorted(wrjobs_content)


# check that colitest run produces correct output
def test_colitest_run(set_vars, run_colitest):
    # check that output/colitest1.cpr is there
    # check that wrdata1/MODEL_STUDY_DONE is there
    colitest_output = set_vars / "output" / "colitest1.cpr"
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    assert colitest_output.is_file()
    assert model_output.is_file()


# compare the output from the colitest job
def test_read_model(set_vars, run_colitest):
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    with open(model_output, "rb") as file:
        data = file.read()
    print("read model file")
    print(data[0:20])
