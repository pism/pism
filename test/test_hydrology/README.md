Run `make exactP` to build the `exactP` extension module.

Apparently the many code warnings generated this way can be ignored.  You will
need the [Cython](http://www.cython.org) library to do this.

Script `simpleP.py` calls the code in `src/verif/tests/exactTestP.c`
through its Python (i.e. Cython) interface:

    $ make exactP
    $ ./simpleP.py

To generate a PISM input file and run test P itself, try

    $ make exactP
    $ ./runTestP.py --pism_path=/path/to/build/directory

To clean up:

    $ make clean

Run `runTestP.py --help` for a summary of its command-line options.

The NetCDF file `inputforP_regression.nc` is used by `test/regression/test_29.py`.

