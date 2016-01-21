Test P
======

Contents of this directory relate to running PISM's Test P which compares
a `-hydrology distributed` run with an exact solution.  This is both a
regression test and a verification test for PISM.

The exact solution is documented in the paper

  * E. Bueler and W. van Pelt (2015) Mass-conserving subglacial hydrology in the Parallel Ice Sheet Model version 0.6. Geoscientific Model Development 8 (6) pp. 1613–1635, [doi:10.5194/gmd-8-1613-2015](http://dx.doi.org/10.5194/gmd-8-1613-2015)

See also the (public) repo containing all materials for that paper, namely [github.com/bueler/hydrolakes](https://github.com/bueler/hydrolakes).


Running the test
----------------

Run `make exactP` to build the `exactP` extension module.

Apparently the many code warnings generated this way can be ignored.  You will
need the [Cython](http://www.cython.org) library to do this.

Scripts `simpleP.py` and `runTestP.py` call the code in
`src/verif/tests/exactTestP.c` through its Python (i.e. Cython) interface:

    $ make exactP
    $ ./simpleP.py

To generate a PISM input file and run test P itself, try

    $ make exactP
    $ ./runTestP.py --pism_path=/path/to/build/directory

To clean up:

    $ make clean

Run `runTestP.py --help` for a summary of its command-line options.

The NetCDF file `inputforP_regression.nc` is used by `test/regression/test_29.py`.

