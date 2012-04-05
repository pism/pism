#include "PISMPNCFile.hh"
#include "PISMNC3File.hh"
#include "PISMNC4File.hh"
#include "PISMProf.hh"
#include "pism_options.hh"

static void compute_Nx_and_Ny(int size, int &Nx, int &Ny) {
  Nx = (int)(0.5 + sqrt(size));

  if (Nx == 0) Nx = 1;

  while (Nx > 0) {
    Ny = size/Nx;
    if (Nx*Ny == size) break;
    Nx--;
  }
}

int main(int argc, char**argv) {

  /* MPI stuff. */
  int mpi_size, mpi_rank, mpi_namelen, ierr;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  string basename = "pism_netcdf_test",
    mode = "netcdf3";

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

  ierr = MPI_Comm_rank(mpi_comm, &mpi_rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(mpi_comm, &mpi_size); CHKERRQ(ierr);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  {
    PISMProf profiler(mpi_comm, mpi_rank, mpi_size);
    int Nx, Ny;

    compute_Nx_and_Ny(mpi_size, Nx, Ny);

    profiler.Nx = Nx;
    profiler.Ny = Ny;

    int event_total = profiler.create("total", "total time"),
      event_alloc = profiler.create("alloc", "time spent allocating memory"),
      event_fill = profiler.create("fill", "time spent filling data"),
      event_output = profiler.create("output", "total output time"),
      event_create = profiler.create("create", "time spent creating the file"),
      event_define = profiler.create("define", "time spent defining variables"),
      event_write = profiler.create("write", "time spent writing data"),
      event_close = profiler.create("close", "time spent closing the file"),
      event_close_reopen = profiler.create("close_reopen", "time spent closing/re-opening the file");

    // Initialize the grid and create non-random data to write:

    bool flag, close_and_reopen;
    // the size of the local block of data:
    int Mx = 801, My = 801, Mz = 201, n_vars = 4;

    ierr = PISMOptionsInt("-Mx", "Number of grid points in the x-direction (for each block)",
                          Mx, flag); CHKERRQ(ierr);
    ierr = PISMOptionsInt("-My", "Number of grid points in the y-direction (for each block)",
                          My, flag); CHKERRQ(ierr);
    ierr = PISMOptionsInt("-Mz", "Number of grid points in the z-direction",
                          Mz, flag); CHKERRQ(ierr);

    ierr = PISMOptionsInt("-n_vars", "Number of variables to write",
                          n_vars, flag); CHKERRQ(ierr);

    ierr = PISMOptionsIsSet("-close", "Close and re-open the file",
                            close_and_reopen); CHKERRQ(ierr);

    if (mpi_rank == 0) {
      printf("each processor allocates ~%lu Mb\n", Mx * My * Mz * sizeof(double) / (1024 * 1024) );
    }

    int start_x = (mpi_rank / Ny),
      start_y = (mpi_rank % Ny);
    // grid parameters:
    double Lx = 1.0 / Nx, Ly = 1.0 / Ny,
      dx = Lx / (Mx - 1), dy = Ly / (My - 1),
      x0 = start_x * Lx, y0 = start_y * Ly,
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

      PISMNCFile *nc = NULL;

      set<string> modes;
      modes.insert("netcdf3"); modes.insert("netcdf4"); modes.insert("pnetcdf");
      ierr = PISMOptionsList(mpi_comm, "-mode", "I/O mode", modes, mode, mode, flag); CHKERRQ(ierr);

      if (mode == "netcdf3") {
        nc = new PISMNC3File(mpi_comm, mpi_rank);
      }
#if (PISM_PARALLEL_NETCDF4==1)
      else if (mode == "netcdf4") {
        nc = new PISMNC4File(mpi_comm, mpi_rank);
      }
#endif
#if (PISM_PNETCDF==1)
      else if (mode == "pnetcdf") {
        nc = new PISMPNCFile(mpi_comm, mpi_rank);
      }
#endif
      else {
        printf("mpi_name: %s rank: %d: mode %s is not supported\n", mpi_name, mpi_rank, mode.c_str());
        PetscFinalize();
        return 1;
      }

      if (nc == NULL) {
        printf("mpi_name: %s rank: %d: memory allocation failed\n", mpi_name, mpi_rank);
        PetscFinalize();
        return 1;
      }

      printf("mpi_name: %s rank: %d: using %s\n", mpi_name, mpi_rank, mode.c_str());

      profiler.begin(event_output);
      {
        profiler.begin(event_create);
        {
          ierr = nc->create(basename + ".nc"); CHKERRQ(ierr);
          int old_fill;
          ierr = nc->set_fill(PISM_NOFILL, old_fill); CHKERRQ(ierr);
        }
        profiler.end(event_create);

        profiler.begin(event_define);
        {
          ierr = nc->def_dim("t", PISM_UNLIMITED); CHKERRQ(ierr);
          ierr = nc->def_dim("x", Mx * Nx); CHKERRQ(ierr);
          ierr = nc->def_dim("y", My * Ny); CHKERRQ(ierr);
          ierr = nc->def_dim("z", Mz); CHKERRQ(ierr);

          vector<string> dims(4);
          dims[0] = "t";
          dims[1] = "x";
          dims[2] = "y";
          dims[3] = "z";

          for (int j = 0; j < n_vars; ++j) {
            char var_name[256];
            snprintf(var_name, 256, "var_%d", j);

            ierr = nc->def_var(var_name, PISM_DOUBLE, dims); CHKERRQ(ierr);
          }

          ierr = nc->enddef(); CHKERRQ(ierr);
        }
        profiler.end(event_define);

        profiler.begin(event_close_reopen);
        if (close_and_reopen) {
          ierr = nc->close(); CHKERRQ(ierr);

          ierr = nc->open(basename + ".nc", PISM_WRITE); CHKERRQ(ierr);
        }
        profiler.end(event_close_reopen);

        profiler.begin(event_write);
        {
          vector<unsigned int> start(4), count(4), imap(4);
          const int T = 0, X = 1, Y = 2, Z = 3;

          start[T] = 0;
          start[X] = Mx * start_x;
          start[Y] = My * start_y;
          start[Z] = 0;

          count[T] = 1;
          count[X] = Mx;
          count[Y] = My;
          count[Z] = Mz;

          for (int j = 0; j < n_vars; ++j) {
            char var_name[256];
            snprintf(var_name, 256, "var_%d", j);

            ierr = nc->put_vara_double(var_name, start, count, data); CHKERRQ(ierr);
          }
        }
        profiler.end(event_write);

        profiler.begin(event_close);
        {
          ierr = nc->close(); CHKERRQ(ierr);
        }
        profiler.end(event_close);
      }
      profiler.end(event_output);

      delete nc;
    }
    profiler.end(event_total);

    ierr = profiler.save_report(basename + "-prof.nc"); CHKERRQ(ierr);

    delete[] data;
  }

  /* Shut down PETSc. */
  PetscFinalize();
  return 0;
}
