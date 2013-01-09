from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os

# If the user set GSL_PREFIX, use it. Otherwise check some standard locations.
prefix = ""
try:
    prefix = os.environ['GSL_PREFIX']
except:
    print "Environment variable GSL_PREFIX not set. Trying known locations..."
    prefixes = ["/usr/", "/usr/local/", "/opt/local/", "/sw/"]

    for path in prefixes:
        print "Checking '%s'..." % path
        try:
            os.stat(path + "include/gsl/gsl_odeiv.h")
            prefix = path
            print "Found GSL in '%s'" % prefix
            break
        except:
            pass

if prefix == "":
    print "Could not find GSL. Stopping..."
    import sys
    sys.exit(1)


# If the user set NO_OPENMP, proceed with these options. Otherwise add options GCC uses.
libraries=['gsl', 'gslcblas']
extra_compile_args=["-O3", "-ffast-math", "-Wall"]

# Define the extension
extension = Extension("exactP",
                      sources=["exactP.pyx",
                               "../../src/verif/tests/exactTestP.c",
                               ],
                      include_dirs=[numpy.get_include(), '../../src/verif/tests/', prefix + '/include'],
                      library_dirs=[prefix + "/lib"],
                      libraries=libraries,
                      extra_compile_args=extra_compile_args,
                      language="c")

setup(
    name = "exactP",
    version = "0.0.1",
    description = "Verification test P helper",
    author = "PISM authors",
    author_email = "help@pism-docs.org",
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extension]
    )
