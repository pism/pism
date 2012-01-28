#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include <netcdf_par.h>

#define CHKERRQ(e) do { \
    if ((e) != NC_NOERR) {                                              \
      printf("Bailing out in file %s, line %d, error:%s.\n",            \
             __FILE__, __LINE__, nc_strerror(e));                       \
      return e;                                                         \
    }                                                                   \
  } while (0)

#define FILE "test_parallel_2d.nc"
#define Y_SIZE 1
#define NDIMS 2

/* This program tries to write data to a 2-dimensional variable using NetCDF's
 * MPI I/O capabilities.
 *
 * Each processor allocates a Y_SIZE * (Y_SIZE + mpi_rank) array
 * and fills it with non-random numbers.
 *
 * "Local" arrays are appended along the "x" dimension.
 *
 * nc_put_varm_double hangs while
 * - writing a two- (or more) dimensional variable
 * - counts along the *first* dimension are different
 * - collective I/O mode is selected
*/

int main(int argc, char **argv)
{
  /* MPI stuff. */
  int mpi_namelen;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  /* Netcdf-4 stuff. */
  int ncid, v1id, dimids[NDIMS];

  size_t start[NDIMS], count[NDIMS];
  ptrdiff_t stride[NDIMS], imap[NDIMS];
  double *data;
  int i, j, stat, local_x_size, local_size, x_size;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  if (argc != 3) {
    if (mpi_rank == 0) {
      printf("Usage: netcdf_parallel_2d vara collective\n"
             "or     netcdf_parallel_2d varm collective\n"
             "or     netcdf_parallel_2d vara independent\n"
             "or     netcdf_parallel_2d varm independent\n");
    }
    /* Shut down MPI. */
    MPI_Finalize();
    return 0;
  }

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  local_x_size = Y_SIZE + mpi_rank;
  local_size = local_x_size * Y_SIZE;

  MPI_Allreduce(&local_x_size, &x_size, 1, MPI_INT, MPI_SUM, comm);

  /* Create a parallel netcdf-4 file. */
  {
    stat = nc_create_par(FILE, NC_NETCDF4|NC_MPIIO, comm, info, &ncid); CHKERRQ(stat);

    /* Create two dimensions. */
    stat = nc_def_dim(ncid, "x", (size_t)x_size, &dimids[0]); CHKERRQ(stat);
    stat = nc_def_dim(ncid, "y", (size_t)Y_SIZE, &dimids[1]); CHKERRQ(stat);

    /* Create one variable. */
    stat = nc_def_var(ncid, "data", NC_DOUBLE, 2, dimids, &v1id); CHKERRQ(stat);

    stat = nc_enddef(ncid); CHKERRQ(stat);
  }

  /* Create phony data. */
  {
    data = (double*)malloc(local_size * sizeof(double));

    if (data == NULL) {
      printf("mpi_rank=%d: Memory allocation failed.\n", mpi_rank);
      exit(1);
    }

    /* Fill with non-random junk to make sure that imap (below) is set correctly. */
    for (i = 0; i < local_x_size; ++i) {
      for (j = 0; j < Y_SIZE; ++j) {
        data[Y_SIZE * i + j] = local_x_size * mpi_rank + i;
      }
    }
  }

  /* Write phony data. */
  {
    start[0] = 0;
    for (j = 0; j < mpi_rank; ++j)
      start[0] += Y_SIZE + j;
    start[1] = 0;

    count[0] = local_x_size;
    count[1] = Y_SIZE;

    /* Note: this mapping from between storage orders in memory and in the file
       is a trivial one. */
    stride[0] = 1;
    stride[1] = 1;

    imap[0] = Y_SIZE;
    imap[1] = 1;

    for (j = 0; j < NDIMS; ++j) {
      printf("mpi_rank=%d: start[%d]=%d count[%d]=%d imap[%d]=%d\n", mpi_rank,
             j, (int)start[j],
             j, (int)count[j],
             j, (int)imap[j]);
    }

    /* Flush buffers (just to make stdout look better). */
    fflush(stdout);
    MPI_Barrier(comm);

    if (strcmp(argv[2], "collective") == 0) {
      stat = nc_var_par_access(ncid, v1id, NC_COLLECTIVE); CHKERRQ(stat);
    } else {
      stat = nc_var_par_access(ncid, v1id, NC_INDEPENDENT); CHKERRQ(stat);
    }

    if (strcmp(argv[1], "varm") == 0) {
      printf("mpi_rank=%d: calling nc_put_varm_double(...)\n", mpi_rank);
      stat = nc_put_varm_double(ncid, v1id, start, count, stride, imap, data); CHKERRQ(stat);
      printf("mpi_rank=%d: exited nc_put_varm_double(...)\n", mpi_rank);
    } else {
      printf("mpi_rank=%d: calling nc_put_vara_double(...)\n", mpi_rank);
      stat = nc_put_vara_double(ncid, v1id, start, count, data); CHKERRQ(stat);
      printf("mpi_rank=%d: exited nc_put_vara_double(...)\n", mpi_rank);
    }

  }

  free(data);

  /* Close the netcdf file. */
  stat = nc_close(ncid); CHKERRQ(stat);

  /* Shut down MPI. */
  MPI_Finalize();
  return 0;
}
