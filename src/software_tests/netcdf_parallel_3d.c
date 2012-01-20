#include <mpi.h>

#include <netcdf_par.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define CHKERRQ(e) do { \
    if ((e) != NC_NOERR) {                                              \
      printf("Bailing out in file %s, line %d, error:%s.\n",            \
             __FILE__, __LINE__, nc_strerror(e));                       \
      return e;                                                         \
    }                                                                   \
  } while (0)

#define FILE "test_parallel_3d.nc"
#define NDIMS 3

int main(int argc, char **argv)
{
  /* MPI stuff. */
  int mpi_namelen;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  /* Netcdf-4 stuff. */
  const size_t BLOCK_LENGTH = 4;
  int ncid, v1id, dimids[NDIMS];

  size_t start[NDIMS], count[NDIMS];
  ptrdiff_t stride[NDIMS], imap[NDIMS];
  double *data;
  int k, stat, local_block_size;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name,
         mpi_size, mpi_rank);

  /* Create a parallel netcdf-4 file. */
  stat = nc_create_par(FILE, NC_NETCDF4|NC_MPIIO, comm,
                       info, &ncid); CHKERRQ(stat);

  /* Create two dimensions. */
  stat = nc_def_dim(ncid, "x", BLOCK_LENGTH, &dimids[0]); CHKERRQ(stat);
  stat = nc_def_dim(ncid, "y", BLOCK_LENGTH, &dimids[1]); CHKERRQ(stat);
  stat = nc_def_dim(ncid, "z", (size_t)mpi_size,     &dimids[2]); CHKERRQ(stat);

  /* Create one var. */
  stat = nc_def_var(ncid, "data", NC_DOUBLE, NDIMS, dimids, &v1id); CHKERRQ(stat);

  stat = nc_enddef(ncid); CHKERRQ(stat);

  /* Set up the slab for this process. */
  start[0] = 0;
  start[1] = 0;
  start[2] = mpi_rank;

  count[0] = BLOCK_LENGTH;
  count[1] = BLOCK_LENGTH;
  count[2] = 1;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  imap[0] = 1;
  imap[1] = BLOCK_LENGTH;
  imap[2] = BLOCK_LENGTH * BLOCK_LENGTH;

  for (k = 0; k < NDIMS; ++k) {
    printf("mpi_rank=%d start[%d]=%d count[%d]=%d imap[%d]=%d\n", mpi_rank,
           k, (int)start[k],
           k, (int)count[k],
           k, (int)imap[k]);
  }

  /* Create phony data. */
  local_block_size = BLOCK_LENGTH * BLOCK_LENGTH;
  data = (double*)malloc(local_block_size * sizeof(double));

  if (data == NULL) {
    printf("mpi_rank=%d: Memory allocation failed.\n", mpi_rank);
    exit(1);
  }

  for (k = 0; k < local_block_size; ++k)
    data[k] = mpi_rank * local_block_size + k;

  stat = nc_var_par_access(ncid, v1id, NC_COLLECTIVE); CHKERRQ(stat);

  /* Write slabs of phony data. */
  stat = nc_put_varm_double(ncid, v1id, start, count, stride, imap, data); CHKERRQ(stat);

  free(data);

  /* Close the netcdf file. */
  stat = nc_close(ncid); CHKERRQ(stat);

  /* Shut down MPI. */
  MPI_Finalize();
  return 0;
}
