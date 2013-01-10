Run `python setup.py build_ext --inplace` to build the `exactP` extension module.
Apparently the many code warnings generated this way can be ignored.  You will
need the Cython library to do this.

Now you can run `simpleP.py` which calls the code in `src/verif/tests/exactTestP.c`
through its Python (i.e. Cython) interface:

    $ python simpleP.py
