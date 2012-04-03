#include "PISMNC4File.hh"
#include "PISMNC3File.hh"

#include <cstdio>		// printf

#define CHKERRQ(e) do { \
    if ((e) != 0) {                                              \
      printf("Bailing out in file %s, line %d.\n", __FILE__, __LINE__); \
      return e;                                                         \
    }                                                                   \
  } while (0)

int main(int argc, char**argv) {

  /* MPI stuff. */
  int mpi_size, mpi_rank, mpi_namelen, ierr;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(mpi_comm, &mpi_size);
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  // Initialize the grid and create non-random data to write:

  // the size of the local block of data:
  int Mx = 11, My = 41, Mz = 11; // these numbers should be all different

  // grid parameters:
  double Lx = 1.0 / mpi_size, Ly = 1.0,
    dx = Lx / (Mx - 1), dy = Ly / (My - 1),
    x0 = dx * Mx * mpi_rank, y0 = 0;

  vector<double> x(Mx), y(My), z(Mz);

  for (int i = 0; i < Mx; ++i)
    x[i] = x0 + dx * i;

  for (int j = 0; j < My; ++j)
    y[j] = y0 + dy * j;

  vector<double> data(Mx * My * Mz);

  for (int i = 0; i < Mx; ++i) {
    for (int j = 0; j < My; ++j) {
      for (int k = 0; k < Mz; ++k) {
        if (k % 2) {
          data[(My * i +  j) * Mz + k] = x[i] + y[j];
        } else {
          data[(My * i +  j) * Mz + k] = mpi_rank;
        }
      }
    }
  }

  // write "data" to a file
  string filename = "pismncfile_test.nc";

  PISMNC4File nc(mpi_comm, mpi_rank);

  ierr = nc.create(filename); CHKERRQ(ierr);

  ierr = nc.def_dim("x", Mx * mpi_size); CHKERRQ(ierr);
  ierr = nc.def_dim("y", My); CHKERRQ(ierr);
  ierr = nc.def_dim("z", Mz); CHKERRQ(ierr);

  vector<string> dims(3);
  dims[0] = "x";
  dims[1] = "y";
  dims[2] = "z";

  ierr = nc.def_var("data", PISM_DOUBLE, dims); CHKERRQ(ierr);

  vector<unsigned int> start(3), count(3), imap(3);

  const int X = 0, Y = 1, Z = 2;

  start[X] = mpi_rank * Mx;
  start[Y] = 0;
  start[Z] = 0;

  count[X] = Mx;
  count[Y] = My;
  count[Z] = Mz;

  imap[X] = My * Mz;
  imap[Y] = Mz;
  imap[Z] = 1;

  ierr = nc.enddef(); CHKERRQ(ierr);

  ierr = nc.put_varm_double("data", start, count, imap, &data[0]); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  /* Shut down MPI. */
  MPI_Finalize();
  return 0;
}
