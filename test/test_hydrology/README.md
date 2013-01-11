Run `make exactP` to build the `exactP` extension module.
Apparently the many code warnings generated this way can be ignored.  You will
need the [Cython](http://www.cython.org) library to do this.

Now you can run `simpleP.py` which calls the code in `src/verif/tests/exactTestP.c`
through its Python (i.e. Cython) interface:

    $ python simpleP.py

To generate a PISM input file and run test P itself, try

    $ ./runTestP.py --pism_path=/path/to/build/directory

Run `runTestP.py --help` for a summary of its command-line options.

The NetCDF file `inputforP_regression.nc` is used by `test/regression/test_29.py`.

