# powr-refactor

This is a development repository for refactoring and improving the numerical stability of the coli part of the [PoWR code](https://github.com/powr-code/PoWR) (Potsdam Wolf-Rayet Stellar Atmospheres). For a description of PoWR and the available models, see [here](https://www.astro.physik.uni-potsdam.de/~wrh/PoWR/powrgrid1.php).

*This library is currently under development!*

## Workflow

The source code is a collection of >400 Fortran77 files with >500 subroutines. For a complete PoWR cycle, different programs are called consecutively. The execution is handled by bash scripts.

In this repository, the development focuses on the coli part of the execution cycle, in particular the [colimo](src/colimo.f) subroutine.

Testing needs to include one and several coli cycles, as well as integration into the PoWR execution cycle.

It also includes compilation of the source code on MacOS and Ubuntu OS with different intel and gnu compilers. The results from different compilers and different optimizations need to be consistent.

## Tasks
- initial set-up of the CI to ensure valid output
- improve the make process: better structure of source directory and make process
- improve the make process: include gnu compiler and debug compiler options
- profile the test runs for time and memory consumption
- implement Fortran best practices in `colimo.f` and modernize to latest standards
- identify computational and numerical bottlenecks and improve the algorithm
