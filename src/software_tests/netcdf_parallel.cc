#include <mpi.h>

extern "C" {
#include <netcdf_par.h>
}

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define CHKERRQ(e) do { \
    if ((e) != NC_NOERR) {                                              \
      printf("Bailing out in file %s, line %d, error:%s.\n", __FILE__, __LINE__, nc_strerror(e)); \
      return e;                                                         \
    }                                                                   \
  } while (0)

#define FILE "test_parallel.nc"
#define NDIMS 2
#define DIMSIZE 12

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
  int data[DIMSIZE*DIMSIZE], i, stat;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  int QTR_DATA = DIMSIZE*DIMSIZE/mpi_size;

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name,
         mpi_size, mpi_rank);

  /* Create a parallel netcdf-4 file. */
  stat = nc_create_par(FILE, NC_NETCDF4|NC_MPIIO, comm,
                       info, &ncid); CHKERRQ(stat);

  /* Create two dimensions. */
  stat = nc_def_dim(ncid, "d1", DIMSIZE, &dimids[0]); CHKERRQ(stat);
  stat = nc_def_dim(ncid, "d2", DIMSIZE, &dimids[1]); CHKERRQ(stat);

  /* Create one var. */
  stat = nc_def_var(ncid, "v1", NC_INT, NDIMS, dimids, &v1id); CHKERRQ(stat);

  stat = nc_enddef(ncid); CHKERRQ(stat);

  /* Set up slab for this process. */
  start[0] = mpi_rank * DIMSIZE/mpi_size;
  start[1] = 0;
  count[0] = DIMSIZE/mpi_size;
  count[1] = DIMSIZE;
  printf("mpi_rank=%d start[0]=%d start[1]=%d count[0]=%d count[1]=%d\n",
         mpi_rank, (int)start[0], (int)start[1], (int)count[0], (int)count[1]);

  /* Create phony data. We’re going to write a 24x24 array of ints,
     in 4 sets of 144. */
  printf("mpi_rank*QTR_DATA=%d (mpi_rank+1)*QTR_DATA-1=%d\n",
         mpi_rank*QTR_DATA, (mpi_rank+1)*QTR_DATA);

  for (i=mpi_rank*QTR_DATA; i<(mpi_rank+1)*QTR_DATA; i++)
    data[i] = mpi_rank;

  // stat = nc_var_par_access(ncid, v1id, NC_COLLECTIVE); CHKERRQ(stat);

  /* Write slabs of phony data. */
  stat = nc_put_vara_int(ncid, v1id, start, count,
                         &data[mpi_rank*QTR_DATA]); CHKERRQ(stat);

  /* Close the netcdf file. */
  stat = nc_close(ncid); CHKERRQ(stat);

  /* Shut down MPI. */
  MPI_Finalize();
  return 0;
}
