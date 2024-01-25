import os
import subprocess
import pytest
import time

import numpy as np


def find_string(str, file, set_vars):
    file = set_vars / "output" / file

    try:
        with open(file, "r") as f:
            output_for_test = f.read()
        return str in output_for_test
    except FileNotFoundError:
        return False


@pytest.fixture()
def run_wrstart(get_chain, set_vars, submit_options):
    if submit_options is None:
        submit_options = ""
    # we actually need to wait for subvarsmit to finish
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

    while True:
        if find_string("WRSTART finished", "wrstart1.log", set_vars) and find_string(
            "wruniq: model finally converged", "wruniq1.log", set_vars
        ):
            break
        time.sleep(1)

    os.system("ls ${POWR_WORK}")
    os.system("ls ${POWR_WORK}/wrdata1")
    os.system("ls ${POWR_WORK}/output")
    os.system("cat ${POWR_WORK}/output/*.cpr")
    yield "ran powr full cycle"
    os.system("rm -rf ${POWR_WORK}/tmp_data")
    return "Cleaned powr tmp data"


def extract_np_between(str, start, end):
    partition = str.partition(start)
    plot = partition[2].partition(end)[0]
    plot_np = np.fromstring(plot, sep=" ")

    return plot_np


@pytest.mark.parametrize("submit_options", ["", " nonopt"])
def test_run_wrstart(
    run_wrstart,
    set_vars,
    get_wrstart1_out_to_match,
    get_wrstart1_cpr_to_match,
    get_wrstart1_plot_to_match,
    get_wruniq1_cpr_to_match,
    get_wruniq1_out_to_match,
    get_wruniq1_plot_to_match,
):
    wrstart1_out_file = set_vars / "output/wrstart1.out"

    with open(wrstart1_out_file, "r") as f:
        output_for_test = f.read()

    for idx, str_searched in enumerate(get_wrstart1_out_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_out_to_match"

    ##########

    wrstart1_cpr_file = set_vars / "output/wrstart1.cpr"

    with open(wrstart1_cpr_file, "r") as f:
        output_for_test = f.read()

    for idx, str_searched in enumerate(get_wrstart1_cpr_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_cpr_to_match"

    # ##########

    wrstart1_plot_file = set_vars / "output/wrstart1.plot"

    with open(wrstart1_plot_file, "r") as f:
        output_for_test = f.read()

    for idx, str_searched in enumerate(get_wrstart1_plot_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wrstart1_plot_to_match"

    #############

    wruniq1_cpr_file = set_vars / "output/wruniq1.cpr"

    with open(wruniq1_cpr_file, "r") as f:
        output_for_test = f.read()

    for idx, str_searched in enumerate(get_wruniq1_cpr_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wruniq1_cpr_to_match"

    ##############

    wruniq1_out_file = set_vars / "output/wruniq1.out"

    with open(wruniq1_out_file, "r") as f:
        output_for_test = f.read()

    for idx, str_searched in enumerate(get_wruniq1_out_to_match):
        assert (
            str_searched in output_for_test
        ), f"String {idx=} not found in get_wruniq1_out_to_match"

    ###############

    wruniq1_plot = set_vars / "output/wruniq1.plot"

    with open(wruniq1_plot, "r") as f:
        plot_for_test = f.read()

    plot_np0 = extract_np_between(plot_for_test, "N=? XYTABLE COLOR=2 PEN=3", "FINISH")
    plot_values0 = np.fromstring(get_wruniq1_plot_to_match[0], sep=" ")
    assert np.allclose(plot_np0, plot_values0, atol=1e-06)

    plot_np1 = extract_np_between(plot_for_test, "N= 1388   PLOTSYMBOL=  5", "FINISH")
    plot_values1 = np.fromstring(get_wruniq1_plot_to_match[1], sep=" ")
    assert np.allclose(plot_np1, plot_values1, atol=1e-06)

    plot_np2 = extract_np_between(plot_for_test, "N=   1388 COLOR=2", "ENDE")
    plot_values2 = np.fromstring(get_wruniq1_plot_to_match[2], sep=" ")
    assert np.allclose(plot_np2, plot_values2, atol=1e-06)

    plot_np3 = extract_np_between(
        plot_for_test,
        "N=     49 COLOR=3 SYMBOL=9 SIZE=0.1",
        "N=     49 COLOR=7 SYMBOL=9 SIZE=0.1",
    )
    plot_values3 = np.fromstring(get_wruniq1_plot_to_match[3], sep=" ")
    assert np.allclose(plot_np3, plot_values3, atol=1e-06)
