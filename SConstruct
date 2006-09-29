import os

petsc_dir = os.environ['PETSC_DIR']
petsc_arch = os.environ['PETSC_ARCH']
mpicc = os.environ['MPICC']
mpicxx = os.environ['MPICXX']
home = os.environ['HOME']
ccflags = "-Wall -Wextra -Wshadow -Wwrite-strings -Wno-unused-parameter"
ccflags += " -Wno-strict-aliasing -Wpointer-arith -Wconversion -Winline"
#ccflags += " -Wcast-qual -Wpadded -Wunreachable-code" # These are excessive
ccflags += " -g3 -DWITH_NETCDF=1 -DWITH_FFTW=1 -DWITH_GSL=1 -pipe"
#ccflags += " -fomit-frame-pointer"
#ccflags += " -pg" # profiling

my_env = Environment(ENV = {'PATH' : os.environ['PATH']},
                     CC=mpicc, CXX=mpicxx,
                     CPPPATH=[os.path.join(petsc_dir, 'include'),
                              os.path.join(petsc_dir, 'bmake', petsc_arch),
                              os.path.join(home, 'usr/include')],
                     LIBPATH=[os.path.join(petsc_dir, 'lib', petsc_arch),
                              '/usr/X11R6/lib', '/usr/lib/atlas/sse2'],
                     RPATH=[os.path.join(petsc_dir, 'lib', petsc_arch)])

deb_env = Environment(ENV = {'PATH' : os.environ['PATH']},
                      CC='mpicc.mpich', CXX='mpicxx.mpich',
                      CPPPATH=['/usr/include/petsc'],
                      LIBPATH=['/usr/X11R6/lib', '/usr/lib/atlas/sse2'],
                      RPATH=[''])

env = my_env
env.Append(CCFLAGS=ccflags)
libice_dir = os.path.join(os.getcwd(), 'obj')
print libice_dir
env['LIBPATH'] += [libice_dir]
env['RPATH'] += [libice_dir]

petsc_libs = ['petsc' + mod for mod in Split('ksp dm mat vec') + [""]]
#petsc_libs += Split('X11 lapack blas mpi_cxx stdc++ dl netcdf_c++ netcdf fftw3 gsl gslcblas')
petsc_libs += Split('X11 lapack blas stdc++ dl netcdf_c++ netcdf fftw3 gsl')

Export('env petsc_libs')

SConscript('src/SConscript', build_dir='obj', duplicate=0)
SConscript('doc/SConscript')

