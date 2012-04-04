#include "PISMPNCFile.hh"
#include "PISMProf.hh"

int main(int argc, char**argv) {

  /* MPI stuff. */
  int mpi_size, mpi_rank, mpi_namelen, ierr;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  string filename = "pism_pnetcdf_test.nc";

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  ierr = MPI_Comm_rank(mpi_comm, &mpi_rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(mpi_comm, &mpi_size); CHKERRQ(ierr);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  {
    PISMProf profiler(mpi_comm, mpi_rank, mpi_size);
    profiler.Nx = mpi_size;
    profiler.Ny = 1;

    int event_total = profiler.create("total", "total time"),
      event_alloc = profiler.create("alloc", "time spent allocating memory"),
      event_fill = profiler.create("fill", "time spent filling data"),
      event_create = profiler.create("create", "time spent creating the file"),
      event_define = profiler.create("define", "time spent defining variables"),
      event_write = profiler.create("write", "time spent writing data"),
      event_close = profiler.create("close", "time spent closing the file");

    // Initialize the grid and create non-random data to write:

    // the size of the local block of data:
    // goal: 1500 x 2800 x 100
    int Mx = 101, My = 401, Mz = 401; // these numbers should be all different

    // grid parameters:
    double Lx = 1.0 / mpi_size, Ly = 1.0,
      dx = Lx / (Mx - 1), dy = Ly / (My - 1),
      x0 = dx * Mx * mpi_rank, y0 = 0,
      *data;

    vector<double> x(Mx), y(My), z(Mz);

    for (int i = 0; i < Mx; ++i)
      x[i] = x0 + dx * i;

    for (int j = 0; j < My; ++j)
      y[j] = y0 + dy * j;

    profiler.begin(event_total);
    {
      profiler.begin(event_alloc);
      {
        data = new double[Mx * My * Mz];

        if (data == NULL) {
          printf("mpi_name: %s rank: %d: memory allocation failed\n", mpi_name, mpi_rank);
          PetscFinalize();
          return 1;
        }
      }
      profiler.end(event_alloc);

      profiler.begin(event_fill);
      {
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
      }
      profiler.end(event_fill);

      // write "data" to a file

      PISMPNCFile nc(mpi_comm, mpi_rank);

      profiler.begin(event_create);
      {
        ierr = nc.create(filename); CHKERRQ(ierr);
        int old_fill;
        ierr = nc.set_fill(PISM_NOFILL, old_fill); CHKERRQ(ierr);
      }
      profiler.end(event_create);

      profiler.begin(event_define);
      {
        ierr = nc.def_dim("x", Mx * mpi_size); CHKERRQ(ierr);
        ierr = nc.def_dim("y", My); CHKERRQ(ierr);
        ierr = nc.def_dim("z", Mz); CHKERRQ(ierr);

        vector<string> dims(3);
        dims[0] = "x";
        dims[1] = "y";
        dims[2] = "z";

        ierr = nc.def_var("data", PISM_DOUBLE, dims); CHKERRQ(ierr);

        ierr = nc.enddef(); CHKERRQ(ierr);
      }
      profiler.end(event_define);

      profiler.begin(event_write);
      {
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

        ierr = nc.put_varm_double("data", start, count, imap, &data[0]); CHKERRQ(ierr);
      }
      profiler.end(event_write);

      profiler.begin(event_close);
      {
        ierr = nc.close(); CHKERRQ(ierr);
      }
      profiler.end(event_close);
    }
    profiler.end(event_total);

    ierr = profiler.save_report(filename + "-prof.nc"); CHKERRQ(ierr);

    delete[] data;
  }

  /* Shut down PETSc. */
  PetscFinalize();
  return 0;
}
