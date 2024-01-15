import os
import subprocess
import pytest
import numpy as np

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


@pytest.fixture()
def run_colitest(get_chain, colitest_options):
    # need aliases and env for this to run on CI
    # run colitest

    # here we need to run with and without nonopt
    if colitest_options is None:
        colitest_options = ""
    colitest_path = "${POWR_WORK}/wrjobs/colitest1"
    colitest_command = colitest_path + colitest_options
    try:
        subprocess.run(
            colitest_command,
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
        "njn_wrh_gen","DEPTH POINTS CONSIDERED:    7  TO   50               CORMAX=66.1896     LOG=  1.82",
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


def extract_np_between(str, start, end):

    partition = str.partition(start)
    plot = partition[2].partition(end)[0]
    plot_np = np.fromstring(plot, sep=' ')

    return plot_np

# check that colitest run produces correct output
@pytest.mark.parametrize("colitest_options", ["", " nonopt"])
def test_colitest(set_vars, get_plot_to_match, run_colitest):

    strs_searched_out = [
        "Maximum Opacity at Depth 1: K= 41042;",
        "Lambda=     303.771;",
        "Opacity= 0.15260E+05;",
        "Opacities set to 0.01 * Background Opacity at    198 frequencies",
        "SMALLPOP=   0.100E-07 chosen different from recommended default   0.100E-11",
        "NEGATIVE BOUND-BOUND COLLISIONAL CROSS SECTION DETECTED (LEVELS: UP=   5, LOW=   2)",
        "MODEL START 21/07/02    16:50:14    70795/0.4D/1600 L=5.3 N=1.5% C=1E-4 Fe=1.4E-3 D4 WNE 10-17",
        "GAMMAC=     80.0   GAMMAL=    800.0",
        "DELTAC=     -1.0   GAMMAR=    800.0   GAMMAD=      0.0",
        "CORRECTIONS REDUCED BY FACTOR 0.50",
    ]

    # check that output/colitest1.cpr is there
    # check that wrdata1/MODEL_STUDY_DONE is there
    colitest_output = set_vars / "output" / "colitest1.cpr"
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    assert colitest_output.is_file()
    assert model_output.is_file()

    # compare the output from the colitest job
    model_output = set_vars / "wrdata1" / "MODEL_STUDY_DONE"
    data_np = np.fromfile(model_output, dtype=float)
    assert len(data_np) == 2246784
    # the original file length is 2233856, so if the coli job does not run
    # the test will fail and the model file length will be 2233856
    assert np.allclose(data_np[256138:256149], MODEL_DATA_REF)

    colitest_file = set_vars / "output/colitest1.out"

    with open(colitest_file, "r") as f:
        output_for_test = f.read()

    for str_searched in strs_searched_out:
        assert str_searched in output_for_test

    ###########################################################################

    colitest_plot = set_vars / "output/colitest1.plot"

    with open(colitest_plot, "r") as f:
        plot_for_test = f.read()

    plot_np0 = extract_np_between(plot_for_test, 'N= 1388   PLOTSYMBOL=  5', 'N=   1388 COLOR=2')
    plot_values0 = np.fromstring(get_plot_to_match[0], sep=' ')
    assert np.allclose(plot_np0, plot_values0, atol=1e-06)

    plot_np1 = extract_np_between(plot_for_test, ' N=   1388 COLOR=2', 'ENDE')
    plot_values1 = np.fromstring(get_plot_to_match[1], sep=' ')
    assert np.allclose(plot_np1, plot_values1, atol=1e-06)

    plot_np2 = extract_np_between(plot_for_test, 'N=     51 COLOR= 2 PEN = 3', 'N=      2 COLOR=3')
    plot_values2 = np.fromstring(get_plot_to_match[2], sep=' ')
    assert np.allclose(plot_np2, plot_values2, atol=1e-06)

    plot_np3 = extract_np_between(plot_for_test, 'N=     50 PEN=4 COLOR=2', 'N=     50 SYMBOL=5 COLOR=4')
    plot_values3 = np.fromstring(get_plot_to_match[3], sep=' ')
    assert np.allclose(plot_np3, plot_values3, atol=1e-06)

    plot_np4 = extract_np_between(plot_for_test, ' N=     50 SYMBOL=5 COLOR=4', 'N=     50 SYMBOL=5 COLOR=9')
    plot_values4 = np.fromstring(get_plot_to_match[4], sep=' ')
    assert np.allclose(plot_np4, plot_values4, atol=1e-06)

    plot_np5 = extract_np_between(plot_for_test, 'N=     50 SYMBOL=5 COLOR=9', 'N=     50 SYMBOL=9 SIZE=0.2')
    plot_values5 = np.fromstring(get_plot_to_match[5], sep=' ')
    assert np.allclose(plot_np5, plot_values5, atol=1e-06)

    plot_np6 = extract_np_between(plot_for_test, 'N=     50 SYMBOL=9 SIZE=0.2', 'N=     50 SYMBOL=10 SIZE=0.2')
    plot_values6 = np.fromstring(get_plot_to_match[6], sep=' ')
    assert np.allclose(plot_np6, plot_values6, atol=1e-06)

    plot_np7 = extract_np_between(plot_for_test, 'N=     50 SYMBOL=10 SIZE=0.2', 'ENDE')
    plot_values7 = np.fromstring(get_plot_to_match[7], sep=' ')
    assert np.allclose(plot_np7, plot_values7, atol=1e-06)

    plot_np12 = extract_np_between(plot_for_test, 'N=   49   PLOTSYMBOL=  1', 'ENDE')
    plot_values12 = np.fromstring(get_plot_to_match[12], sep=' ')
    assert np.allclose(plot_np12, plot_values12, atol=1e-06)
