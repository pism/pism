#!/usr/bin/env python

from sys import argv, stderr, path, exit
from os import getenv, putenv, system

# get the path to PISM's Python bindings
PISM_PATH = argv[1]

# get the path to Python scripts to run
script_path = argv[2]

# get the test number
N = int(argv[3])

# set the PYTHONPATH environment variable so that Python scripts run with
# os.system can find PISM's Python bindings
putenv("PYTHONPATH", "%s:%s" % (PISM_PATH, getenv("PYTHONPATH")))

def run_test(name):
    stderr.write("Running %s\n" % name)
    err = system("%s/%s" % (script_path, name))
    if err != 0:
        exit(-1)

# Check if we can 'import PISM'
try:
    path.append(PISM_PATH)
    import PISM
except:
    stderr.write("ERROR: can't import PISM. Did you run 'make install'?\n")
    exit(-1)

# Run a test
if   N == 1:
    run_test("ssa_tests/ssa_testi.py -Mx 3 -My 61")
elif N == 2:
    run_test("ssa_tests/ssa_testj.py")
elif N == 3:
    run_test("ssa_tests/ssa_test_linear.py")
elif N == 4:
    run_test("ssa_tests/ssa_test_plug.py -Mx 3 -My 61")
