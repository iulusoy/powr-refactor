import os
import subprocess
import pytest
import numpy as np

import re

MODEL_DATA_REF = np.array(
    [
        200.48528618028934,
        170.96813898862962,
        141.27907698674588,
        113.21541313779215,
        88.20625902763122,
        67.06805146763016,
        49.991012566374565,
        36.694735175459314,
        26.638377911447826,
        19.19939157825736,
        13.787069162858446,
    ]
)


@pytest.fixture(scope="module")
def run_colitest(get_chain):
    # need aliases and env for this to run on CI
    # run colitest
    colitest_command = "${POWR_WORK}/wrjobs/colitest1"

    try:
        temp = subprocess.run(
            colitest_command,
            shell=True,
            check=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
        )
    except Exception as error:
        print(temp.stderr)
        print(error)

    print(temp.stdout)

    os.system("cat ${POWR_WORK}/output/colitest1.cpr")
    yield "ran colitest"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    os.system("rm -rf ${POWR_WORK}/tmp_2day/*")
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
        "FEDAT",
        "FEDAT_FORMAL",
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
    data_np = np.fromfile(model_output, dtype=float)
    assert len(data_np) == 2246784
    assert np.allclose(data_np[256138:256149], MODEL_DATA_REF)


    colitest_file_reference = set_vars / "../test/data/colitest1.out"
    colitest_file = set_vars / "output/colitest1.out"

    with open(colitest_file_reference, "r") as f:
        output = f.read()

    with open(colitest_file, "r") as f:
        output_for_test = f.read()

    out1 = output.find("Maximum Opacity at Depth 1: K= 41042;")
    out2 = output_for_test.find("Maximum Opacity at Depth 1: K= 41042;")
    assert out1 == out2

    out1 = output.find("Lambda=     303.771;")
    out2 = output_for_test.find("Lambda=     303.771;")
    assert out1 == out2

    out1 = output.find("Opacity= 0.15260E+05;")
    out2 = output_for_test.find("Opacity= 0.15260E+05;")
    assert out1 == out2

    out1 = output.find("Opacities set to 0.01 * Background Opacity at    198 frequencies")
    out2 = output_for_test.find("Opacities set to 0.01 * Background Opacity at    198 frequencies")
    assert out1 == out2

    out1 = output.find("SMALLPOP=   0.100E-07 chosen different from recommended default   0.100E-11")
    out2 = output_for_test.find("SMALLPOP=   0.100E-07 chosen different from recommended default   0.100E-11")
    assert out1 == out2

    out1 = output.find("NEGATIVE BOUND-BOUND COLLISIONAL CROSS SECTION DETECTED (LEVELS: UP=   5, LOW=   2)")
    out2 = output_for_test.find("NEGATIVE BOUND-BOUND COLLISIONAL CROSS SECTION DETECTED (LEVELS: UP=   5, LOW=   2)")
    assert out1 == out2

    out1 = output.find("MODEL START 21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17")
    out2 = output_for_test.find("MODEL START 21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17")
    assert out1 == out2

    out1 = output.find("MAX: 67.1896  (G  4....1e  L=  7)       2ND: 65.5906  (G  4....1e  L=  8)")
    out2 = output_for_test.find("MAX: 67.1896  (G  4....1e  L=  7)       2ND: 65.5906  (G  4....1e  L=  8)")
    assert out1 == out2

    out1 = output.find("MIN:  0.5000  (G  5....6e  L=  7)       2ND:  0.5000  (G  5....6e  L=  8)")
    out2 = output_for_test.find("MIN:  0.5000  (G  5....6e  L=  7)       2ND:  0.5000  (G  5....6e  L=  8)")
    assert out1 == out2

    out1 = output.find("DEPTH POINTS CONSIDERED:    7  TO   50               CORMAX=66.1896     LOG=  1.82")
    out2 = output_for_test.find("DEPTH POINTS CONSIDERED:    7  TO   50               CORMAX=66.1896     LOG=  1.82")
    assert out1 == out2

    out1 = output.find("GAMMAC=     80.0   GAMMAL=    800.0")
    out2 = output_for_test.find("GAMMAC=     80.0   GAMMAL=    800.0")
    assert out1 == out2

    out1 = output.find("DELTAC=     -1.0   GAMMAR=    800.0   GAMMAD=      0.0")
    out2 = output_for_test.find("DELTAC=     -1.0   GAMMAR=    800.0   GAMMAD=      0.0")
    assert out1 == out2

    out1 = output.find("CORRECTIONS REDUCED BY FACTOR 0.50")
    out2 = output_for_test.find("CORRECTIONS REDUCED BY FACTOR 0.50")
    assert out1 == out2

    #######################################################################

    colitest_plot_reference = set_vars / "../test/data/colitest1.plot"
    colitest_plot = set_vars / "output/colitest1.plot"

    with open(colitest_plot_reference, "r") as f:
        output = f.read()

    with open(colitest_plot, "r") as f:
        output_for_test = f.read()

    out1 = output.find("0.69897000      0.72478278      0.74914704      0.77221675      0.79412257")
    out2 = output_for_test.find("0.69897000      0.72478278      0.74914704      0.77221675      0.79412257")
    assert out1 == out2

    out1 = output.find("-21.832692      -22.997487      -24.157987")
    out2 = output_for_test.find("-21.832692      -22.997487      -24.157987")
    assert out1 == out2

    out1 = output.find("0.69897000      0.72478278      0.74914704      0.77221675      0.79412257 ")
    out2 = output_for_test.find("0.69897000      0.72478278      0.74914704      0.77221675      0.79412257 ")
    assert out1 == out2

    out1 = output.find("PLOT: HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    out2 = output_for_test.find("PLOT: HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    assert out1 == out2

    out1 = output.find("PLOT   :HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    out2 = output_for_test.find("PLOT   :HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    assert out1 == out2

    out1 = output.find("PLOT   :HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    out2 = output_for_test.find("PLOT   :HSUM: M21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17 JOB No.     54")
    assert out1 == out2
