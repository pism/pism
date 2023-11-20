To build the PISM python bindings, you need to have both SWIG and petsc4py installed.

To build the bindings, you will need to first set the cmake flag Pism_BUILD_PYTHON_BINDINGS.  In addition,
the build process requires paths to various pieces of python as well to petsc4py.
We try to detect these pieces automatically, but this may fail and you may need to set the
following advanced cmake variables explicitly (with examples in parentheses).

PYTHON_EXECUTABLE (/opt/local/bin/python)
PYTHON_INCLUDES (/opt/local/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6)
PYTHON_LIBRARY  (/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/config/libpython2.6.dylib)
PETSC4PY_INCLUDES (/Users/david/.local/lib/python2.6/site-packages/petsc4py/include)

Try a 'make' and see if it builds.  There will be a flurry of warnings; go ahead and ignore them.

After having done a 'make install', the python libraries will be found in something like

${CMAKE_INSTALL_PREFIX}/lib/python2.6/site-packages/PISM/...etc...

You will need to add this location to your PYTHONPATH, e.g.

export PYTHONPATH=${PYTHONPATH}:/usr/local/lib/python2.6/site-packages

NOTE: The path terminates with site-packages, not PISM.

Now give it a shot.  Simply try "python -c 'import PISM'" to see if you get an error message.

If this works, this is very promising.  Now go to the examples/python/ssa_tests/
directory.  Try

  python ssa_testj.py

  python ssa_testj.py -Mx 121 -My 121 -snes_monitor

  mpiexec -np 4 python ssa_testj.py
