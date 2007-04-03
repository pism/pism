import os

petsc_dir = os.environ.get('PETSC_DIR', '/usr/lib/petsc')
petsc_arch = os.environ.get('PETSC_ARCH', 'linux-gnu-c-opt')
mpicc = os.environ.get('MPICC', 'mpicc')
mpicxx = os.environ.get('MPICXX', 'mpicxx')
home = os.environ.get('HOME')
ccflags = "-Wall -Wextra -Wshadow -Wwrite-strings -Wno-unused-parameter"
ccflags += " -Wno-strict-aliasing -Wpointer-arith -Wconversion -Winline"
#ccflags += " -Wcast-qual -Wpadded -Wunreachable-code" # These are excessive
ccflags += " -g3 -DWITH_FFTW=1 -DWITH_GSL=1 -pipe"
#ccflags += " -fomit-frame-pointer"
#ccflags += " -pg" # profiling

my_env = Environment(ENV = {'PATH' : os.environ['PATH']},
                     CC=mpicc, CXX=mpicxx,
                     CPPPATH=[os.path.join(petsc_dir, 'include'),
                              os.path.join(petsc_dir, 'bmake', petsc_arch)],
                     LIBPATH=[os.path.join(petsc_dir, 'lib', petsc_arch),
                              '/usr/X11R6/lib', '/usr/lib/atlas/sse2'],
                     RPATH=[''])

deb_env = Environment(ENV = {'PATH' : os.environ['PATH']},
                      CC='mpicc.mpich', CXX='mpicxx.mpich',
                      CPPPATH=['/usr/include/petsc'],
                      LIBPATH=['/usr/X11R6/lib', '/usr/lib/atlas/sse2'],
                      RPATH=[''])

debug = ARGUMENTS.get('debug', 0)
if int(debug):
    deb_env.Append(LIBPATH = ['/usr/lib/petsc/lib/debug'])

env = my_env

conf = Configure(env)

if conf.CheckCHeader('netcdf.h'):
    ccflags += ' -DWITH_NETCDF=1'
else:
    print 'No netCDF support.'
    ccflags += ' -DWITH_NETCDF=0'

env = conf.Finish()

# Get the necessary libraries and linker flags from the petscconf file.
pcc_linker_libs = filter(lambda s: s[:15] == 'PCC_LINKER_LIBS',
                          open(os.path.join(petsc_dir, 'bmake', petsc_arch, 'petscconf')))[0].split(' ')
                          #open('/home/jed/nelchina/petscconf_NELCHINA'))[0].split(' ')
rpath = filter(lambda s: s[:11] == '-Wl,-rpath,', pcc_linker_libs)
libpath = filter(lambda s: s[:2] == '-L', pcc_linker_libs)

print 'rpath = ' + ' '.join(rpath)
print 'libpath = ' + ' '.join(libpath)
if (0):
    env.Append(RPATH = [s[11:] for s in rpath])
    env.Append(LIBPATH = [s[2:] for s in libpath])
else:
    print 'not using rpath or libpath'

env.Append(CCFLAGS=ccflags)
libpism_dir = os.path.join(os.getcwd(), 'obj')
print libpism_dir
env.Append(LIBPATH = [libpism_dir])
env.Append(RPATH = [libpism_dir])

petsc_libs = ['petsc' + mod for mod in Split('ksp dm mat vec') + [""]]
pism_libs = petsc_libs + Split('X11 lapack blas stdc++ dl netcdf fftw3 gsl')

Export('env pism_libs')

SConscript('src/SConscript', build_dir='obj', duplicate=0)
SConscript('doc/SConscript', duplicate=0)

