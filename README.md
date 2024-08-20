# powr-refactor

![License: GPL-3](https://img.shields.io/github/license/ssciwr/powr-refactor)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/ssciwr/powr-refactor/ci.yml?branch=main)
![Language](https://img.shields.io/github/languages/top/ssciwr/powr-refactor)

This is a development repository for refactoring and improving the numerical stability of the coli part of the [PoWR code](https://github.com/powr-code/PoWR) (Potsdam Wolf-Rayet Stellar Atmospheres). For a description of PoWR and the available models, see [here](https://www.astro.physik.uni-potsdam.de/~wrh/PoWR/powrgrid1.php).

*This library is currently under development!*

## Workflow

The source code is a collection of >400 Fortran77 files with >500 subroutines. For a complete PoWR cycle, different programs are called consecutively. The execution is handled by bash scripts that call each other. These bash scripts are placed in `powr/dummychain` and `powr/proc.dir`.

In this repository, the development focuses on the coli part of the execution cycle, in particular the [colimo](src/colimo.f) subroutine.

Testing needs to include one and several coli cycles, as well as integration into the PoWR execution cycle. For this, we run `colitest` and `wrstart`. The integration tests are only run before merging into the main branch, as they are quite costly in terms of compute time.

Testing also also includes compilation of the source code on MacOS and Ubuntu OS with different intel and gnu compilers. The results from different compilers and different optimizations need to be consistent.

## Tasks
- [x] initial set-up of the CI to ensure valid output
- [x] improve the make process: better structure of source directory and make process
- [x] improve the make process: include ifx and gnu compiler and debug compiler options
- profile the test runs for time and memory consumption
- implement Fortran best practices in `colimo.f` and modernize to latest standards
- identify computational and numerical bottlenecks and improve the algorithm

# Compiling PoWR

To compile PoWR, make sure you have a Fortran compiler and math libraries installed. We recommend intel Fortran and MKL libraries. These can be obtained through the [OneAPI toolkits](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html).
After installing the compiler(s) and libraries, make sure that the environment variables in your shell are set correctly - so that you can compile, link and execute the code. Intel OneAPI provides a `set_up_vars.sh` script that should (in principle, sometimes there can be version conflicts if you have several sets of compilers installed) set up everything for you.

To compile, simply execute `make` in the base repository directory of PoWR. This will provide you with (optimized) executables in the `powr/exe.dir` directory, compiled with `ifx`, the modern Fortran compiler by Intel.

Other compilation options are:
- `make small`: Compilation with `ifx` and optimization, of only `coli` and `steal` as targets
- `make debug`: Compilation with debug options and no optimization for `ifx`, only `coli` and `steal` are compiled and the executables are placed in `powr/exe_dev.dir`
- `make debug_all`: Compilation with debug options and no optimization for `ifx`, all programs are compiled and the executables are placed in `powr/exe_dev.dir`
- `make intel_classic`: Compilation with the classic `ifort` compiler with optimization, all programs are compiled and the executables are placed in `powr/exe.dir`
- `make intel_classic_debug`: Compilation with debug options and no optimization with the classic `ifort` compiler, only `coli` and `steal` are compiled and the executables are placed in `powr/exe_dev.dir`
- `make gfortran`: Compilation with the gnu Fortran compiler, only `coli` is compiled
- `make clean`: remove all binaries, object and module files
- `make clean_build`: remove all object files to allow recompilation with/without debug options, without removing the compiled binaries in the `exe.dir`/`exe_dev.dir` folders

The object and module files are placed in `build`, module files in `modules`, library files in `lib`, to allow a clean separation of source files, object and module files, and binary executables.

*Note: In the future, the make process should make heavier use of encapsulated module files.*

# Tests

## Detailed description of the tests

The tests are driven by pytest. The configuration of the tests exports all the necessary variables and sets the stage for the bash scripts, that are then called in subprocesses. The output of the different jobs is compared to reference output in assert statements.

There are also testmodel files included in the repository. These include runs that are numerically unstable and need to be debugged. The testmodel files can be downloaded using
```
wget -O testmodels.tgz https://heibox.uni-heidelberg.de/f/a62c7ae5559d43a0a8b2/?dl=1
```

### `coli` test
The coli test runs follow this workflow:
- Create a new chain, ie 1, by `makechain 1`. This copies some scripts and executables into different folders:

| Folder      | Purpose |
| ----------- | ----------- |
| wrdata1     | data directory of the current results |
| scratch subdirectories | directories to handle execution cycle process |
| wrjobs | collection of scripts to run |
| output | directory containing the global output |
| tmp_data | intermediate results? commonly scratch? |
| tmp_2day | ? |

The most important files are in the `wrdata1` directory: `CARDS, DATOM, FEDAT, FEDAT_FORMAL, FGRID, FORMAL_CARDS, MODEL, MODEL_STUDY, NEWDATOM_INPUT, NEWFORMAL_CARDS_INPUT, next_job, next_jobz`  
Most of these are input files, some control the job execution.

- The test is then run by `sub colitest1` (or directly by calling colitest1). This creates the run log in `output` (`colitest1.log` and `colitest1.cpr`) - these are checked if the run are successful (`COLITEST  finished` in `log`). The results in `cpr` are compared to the reference file.

The output that is generated in `wrdata1` is also compared: `MODEL_STUDY_DONE, MODEL_STUDY_DONE_STEAL`


### Integration test

- First, a new chain is generated using `makechain 1`.
- The integration test calls `wrstart` through `submit`. `wrstart` then calls `wruniq` which handles the COMO / COLI / STEAL program cycles. Since all of these processes are detached from the initial submit process, we need to check for completion by regularly parsing the log files.

### Unit tests
TBD

### Testing through Python
`colitest` (small integration test, coli and steal) and the full integration test (including the full cycle of program executions) are run using `pytest`. The Python test session set-up and tear-down are included in `tests`.

### FAQ
- `ifx` is returning a runtime error when executing coli: This may happen if the environment variables for oneapi are not set correctly. Force-source the set-up using
```
source /opt/intel/oneapi/2024.2/oneapi-vars.sh --force
```