import os
import tempfile
import string

prefix = ARGUMENTS.get('prefix', os.getcwd());
petsc_dir = os.environ.get('PETSC_DIR', '/usr/lib/petsc')
petsc_arch = os.environ.get('PETSC_ARCH', 'linux-gnu-c-opt')
mpicc = os.environ.get('MPICC', 'mpicc')
mpicxx = os.environ.get('MPICXX', 'mpicxx')
home = os.environ.get('HOME')
my_ccflags = "-Wall -Wextra -Wshadow -Wwrite-strings -Wno-unused-parameter"
my_ccflags += " -Wno-strict-aliasing -Wpointer-arith -Wconversion -Winline"
#my_ccflags += " -Wcast-qual -Wpadded -Wunreachable-code" # These are excessive
my_ccflags += " -g3 -pipe"

mt = os.environ.get('MARGIN_TRICK', '0')
mt2 = os.environ.get('MARGIN_TRICK_TWO', '0')
defines = ' -DMARGIN_TRICK=' + mt + ' -DMARGIN_TRICK_TWO=' + mt2

makefile_contents = ('include ' + os.path.join(petsc_dir, 'bmake/common/base') + '\n'
            + 'ALL : \n'
            + '\t@echo ${PETSC_COMPILE_SINGLE}\n'
            + '\t@echo ${CLINKER}\n'
            + '\t@echo ${PETSC_LIB}')

#print makefile_contents

makefile = tempfile.mkstemp()
os.write(makefile[0], makefile_contents)
outputfile = tempfile.mkstemp()
os.system('make ALL -f ' + makefile[1] + ' > ' + outputfile[1])
os.remove(makefile[1])
output = open(outputfile[1])
petsc_compile_single = output.readline()
petsc_clinker = output.readline()
petsc_lib = output.readline()[:-1]
os.remove(outputfile[1])

#print petsc_compile_single, petsc_clinker, petsc_lib
def xworks(string):
    return os.system(string + ' 2>\&1 > /dev/null') == 0
petsc_cc = petsc_compile_single.split(' ')[0]
petsc_cxx = petsc_cc.replace('mpicc', 'mpicxx')
if xworks(petsc_cxx + ' -v'):
    print 'Using mpicxx: ', petsc_cxx
else:
    print 'Cannot use mpicxx: ', petsc_cxx
    petsc_cxx = petsc_cc.replace('mpicc', 'mpic++')
    if xworks(petsc_cxx + ' -v'):
        print 'Using mpic++: ', petsc_cxx
    else:
        print 'Cannot use mpic++: ', petsc_cxx
        print 'Will try with mpicc: ', petsc_cc
        petsc_cxx = petsc_cc

petsc_ccflags = [w for w in petsc_compile_single.split(' ')[2:]
                if w[:2] in ('-f', '-W', '-g', '-p')]
petsc_cpppath = [w[2:] for w in petsc_compile_single.split(' ')[2:]
                 if w[:2] == '-I']
petsc_libpath = [w[2:] for w in petsc_lib.split(' ')
                 if w[:2] == '-L']
petsc_rpath = [w[11:] for w in petsc_lib.split(' ')
               if w[:11] == '-Wl,-rpath,']
petsc_libs = [w[2:] for w in petsc_lib.split(' ')
              if w[:2] == '-l']

petsc_env = Environment(ENV = {'PATH' : os.environ['PATH']},
                        CC=petsc_cc, CXX=petsc_cxx,
                        CPPPATH=petsc_cpppath,
                        LIBPATH=petsc_libpath,
                        RPATH=petsc_rpath)

env = petsc_env
ccflags = string.join(petsc_ccflags)

conf = Configure(env)
if not conf.CheckLibWithHeader('netcdf', 'netcdf.h', 'C'):
    print 'Pism needs the netCDF (>=3.6.0) library and header.'
    Exit(1)
if not conf.CheckLibWithHeader(Split('gsl gslcblas'), 'gsl/gsl_integration.h', 'C'):
    print 'Pism needs the GSL library and headers.'
    Exit(1)
if conf.CheckLibWithHeader('fftw3', 'fftw3.h', 'C'):
    defines += " -DWITH_FFTW=1"
else:
    defines += " -DWITH_FFTW=0"
env = conf.Finish()

libpism_dir = os.path.join(prefix, 'lib')
#print libpism_dir
env.Append(CCFLAGS = ccflags + defines)
env.Append(LIBPATH = ['.', libpism_dir])
env.Append(RPATH = ['.', libpism_dir])

#petsc_libs = ['petsc' + mod for mod in Split('ksp dm mat vec') + [""]]
pism_libs = petsc_libs + Split('stdc++ netcdf fftw3 gsl gslcblas')

if False:
    print 'CC        = ', env['CC']
    print 'CXX       = ', env['CXX']
    print 'CCFLAGS   = ', env['CCFLAGS']
    print 'CPPPATH   = ', env['CPPPATH']
    print 'LIBPATH   = ', env['LIBPATH']
    print 'RPATH     = ', env['RPATH']
    print 'pism_libs = ', pism_libs

Export('env prefix pism_libs')

SConscript('src/SConscript')
SConscript('doc/SConscript')

