import os
import subprocess
import pytest
import time
import re

import numpy as np


def find_string(str, file, set_vars):
    file = set_vars / "output" / file

    try:
        with open(file, "r") as f:
            output_for_test = f.read()
        return str in output_for_test
    except FileNotFoundError:
        return False


def remove_string_get_array(mystring):
    # remove all strings from the list of numbers
    numbers = re.sub("[^0-9. \-]", " ", mystring)  # noqa
    myarray = np.fromstring(numbers, sep=" ")
    return myarray


def extract_np_between(str, start, end):
    partition = str.partition(start)
    mystring = partition[2].partition(end)[0]
    myarray = remove_string_get_array(mystring)
    return myarray


@pytest.fixture(scope="module")
def run_wrstart(get_chain, set_vars):
    # we actually need to wait for subvarsmit to finish
    # wrstart finishes first then wruniq
    # for wruniq, in wruniq1.log there appears a line
    # ""
    print("Starting wrstart run")
    submit_path = "${POWR_WORK}/proc.dir/submit.com wrstart1"
    submit_command = submit_path
    # run wrstart

    try:
        result = subprocess.run(
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
    print("wrstart stdout")
    print(result.stdout)
    print("wrstart stderr")
    print(result.stderr)
    print("done with wrstart run")

    while True:
        if find_string("WRSTART finished", "wrstart1.log", set_vars) and find_string(
            "wruniq: model finally converged", "wruniq1.log", set_vars
        ):
            break
        time.sleep(1)
        print("Looking for finished job..")
        print(result.stdout)
        print(result.stderr)
        os.system(f"ls {set_vars}/output")
        print("cat wrstart1.log")
        os.system(f"cat {set_vars}/output/wrstart1.log")
        print("cat wruniq1.log")
        os.system(f"cat {set_vars}/output/wruniq1.log")

    os.system("ls ${POWR_WORK}")
    os.system("ls ${POWR_WORK}/wrdata1")
    os.system("ls ${POWR_WORK}/output")
    os.system("cat ${POWR_WORK}/output/*.cpr")
    yield "ran powr full cycle"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    return "Cleaned powr tmp data"


def test_wrstart_out(run_wrstart, set_vars, get_wrstart1_out_to_match):
    wrstart1_out_file = set_vars / "output/wrstart1.out"
    with open(wrstart1_out_file, "r") as f:
        output_for_test = f.read()
    for idx, str_searched in enumerate(get_wrstart1_out_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_out_to_match"


def test_wrstart_cpr(run_wrstart, set_vars, get_wrstart1_cpr_to_match):
    wrstart1_cpr_file = set_vars / "output/wrstart1.cpr"
    with open(wrstart1_cpr_file, "r") as f:
        output_for_test = f.read()
    for idx, str_searched in enumerate(get_wrstart1_cpr_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_cpr_to_match"


def test_wrstart_plot(run_wrstart, set_vars, get_wrstart1_plot_to_match):
    wrstart1_plot_file = set_vars / "output/wrstart1.plot"
    with open(wrstart1_plot_file, "r") as f:
        output_for_test = f.read()
    for idx, str_searched in enumerate(get_wrstart1_plot_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_plot_to_match"


def test_wruniq_cpr(run_wrstart, set_vars, get_wruniq1_cpr_to_match):
    wruniq1_cpr_file = set_vars / "output/wruniq1.cpr"
    with open(wruniq1_cpr_file, "r") as f:
        output_for_test = f.read()
    for idx, str_searched in enumerate(get_wruniq1_cpr_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wruniq1_cpr_to_match"


def test_wruniq_out(run_wrstart, set_vars, get_wruniq1_out_to_match):
    wruniq1_out_file = set_vars / "output/wruniq1.out"
    with open(wruniq1_out_file, "r") as f:
        output_for_test = f.read()
    for idx, str_searched in enumerate(get_wruniq1_out_to_match[0]):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wruniq1_out_to_match"

    # test for float values
    wruniq1_out_np0 = extract_np_between(
        output_for_test,
        "(ATOMS PER CM+3)    (EL. PER CM+3)    (NLTE)       (NLTE)",
        "RADII AND TEMPERATURES FOR DIFFERENT OPTICAL DEPTHS",
    )
    out_values0 = remove_string_get_array(get_wruniq1_out_to_match[1][0])
    assert np.allclose(wruniq1_out_np0, out_values0, atol=1e-06)

    wruniq1_out_np1 = extract_np_between(
        output_for_test,
        "(KM/S)           (KM/S PER RSTAR)  DENSCON     FILLFAC    MASS (XMU)",
        " ------------------------------------ 1          MODEL START",
    )
    out_values1 = remove_string_get_array(get_wruniq1_out_to_match[1][1])
    assert np.allclose(wruniq1_out_np1, out_values1, atol=1e-06)

    wruniq1_out_np2 = extract_np_between(
        output_for_test,
        "(KELVIN)    (KELVIN)    B(TEFF)    (ERG/CM2/SEC)        (MAGNITUDES)",
        "EFFECTIVE TEMPERATURE:  70581. KELVIN",
    )
    out_values2 = remove_string_get_array(get_wruniq1_out_to_match[1][2])
    # don't check the large values in the end
    wruniq1_out_np2 = wruniq1_out_np2[:2282]
    out_values2 = out_values2[:2282]
    assert np.allclose(wruniq1_out_np2, out_values2, atol=1e-6)


def test_wruniq_plot(run_wrstart, set_vars, get_wruniq1_plot_to_match):
    wruniq1_plot = set_vars / "output/wruniq1.plot"
    with open(wruniq1_plot, "r") as f:
        plot_for_test = f.read()

    plot_np0 = extract_np_between(plot_for_test, "N=? XYTABLE COLOR=2 PEN=3", "FINISH")
    plot_values0 = remove_string_get_array(get_wruniq1_plot_to_match[0])
    assert np.allclose(plot_np0, plot_values0, atol=1e-06)

    plot_np1 = extract_np_between(
        plot_for_test, "N= 1388   PLOTSYMBOL=  5", "N=   1388 COLOR=2"
    )
    plot_values1 = remove_string_get_array(get_wruniq1_plot_to_match[1])
    assert np.allclose(plot_np1, plot_values1, atol=1e-06)

    plot_np2 = extract_np_between(plot_for_test, "N=   1388 COLOR=2", "ENDE")
    plot_values2 = remove_string_get_array(get_wruniq1_plot_to_match[2])
    assert np.allclose(plot_np2, plot_values2, atol=1e-06)

    plot_np3 = extract_np_between(
        plot_for_test,
        "N=     49 COLOR=3 SYMBOL=9 SIZE=0.1",
        "N=     49 COLOR=7 SYMBOL=9 SIZE=0.1",
    )
    plot_values3 = remove_string_get_array(get_wruniq1_plot_to_match[3])
    assert np.allclose(plot_np3, plot_values3, atol=1e-06)
