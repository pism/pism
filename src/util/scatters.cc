/*!
 * This code sets up scatters of a distributed Vec to a given rank to process data using
 * serial code.
 *
 * It is *not* used at the moment.
 */
#if 0
struct VecAndScatter {
  VecScatter scatter;
  Vec v;
};

/*!
 * Allocate the scatter from a part of a parallel Vec to a target rank.
 *
 * The caller is responsible for de-allocating both the scatter and the target Vec.
 */
static VecAndScatter scatter_part(Vec v_in, int start, int length, int target_rank) {
  PetscErrorCode ierr;
  int rank;
  VecAndScatter result;
  IS is;

  MPI_Comm_rank(PetscObjectComm((PetscObject)v_in), &rank);

  if (rank != target_rank) {
    length = 0;
  }

  ierr = VecCreateSeq(PETSC_COMM_SELF, length, &result.v);
  PISM_CHK(ierr, "VecCreateSeq");

  ierr = ISCreateStride(PETSC_COMM_SELF, length, start, 1, &is);
  PISM_CHK(ierr, "ISCreateStride");

  ierr = VecScatterCreate(v_in, is, result.v, NULL, &result.scatter);
  PISM_CHK(ierr, "VecScatterCreate");

  ierr = ISDestroy(&is);
  PISM_CHK(ierr, "ISDestroy");

  return result;
}

/*!
 * Allocate a natural Vec for a given DM.
 *
 * The caller is responsible for de-allocating the Vec returned by this function.
 */
static Vec get_natural_work(DM dm) {
  PetscErrorCode ierr;
  Vec result;

  ierr = PetscObjectQuery((PetscObject)dm, "natural_work", (PetscObject*)&result);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (result == NULL) {
    Vec v = NULL;
    ierr = DMDACreateNaturalVector(dm, &v);
    PISM_CHK(ierr, "DMDACreateNaturalVector");

    ierr = PetscObjectCompose((PetscObject)dm, "natural_work", (PetscObject)(v));
    PISM_CHK(ierr, "PetscObjectCompose");

    result = v;

    ierr = VecDestroy(&v);
    PISM_CHK(ierr, "VecDestroy");
  }

  return result;
}

/*!
 * Given a DM, allocate a rank 0 target Vec that can be used to gather a part of a
 * "global" Vec on rank 0. Arguments "start" and "length" define the part in question.
 *
 * The caller is responsible for de-allocating the Vec returned by this function.
 */
static Vec proc0_copy(DM dm, int start, int length) {
  Vec v_proc0 = NULL;
  PetscErrorCode ierr = 0;

  ierr = PetscObjectQuery((PetscObject)dm, "v_proc0", (PetscObject*)&v_proc0);
  PISM_CHK(ierr, "PetscObjectQuery");
                                                                                          ;
  if (v_proc0 == NULL) {

    // natural_work will be destroyed at the end of scope, but it will
    // only decrement the reference counter incremented by
    // PetscObjectCompose below.
    auto natural_work = get_natural_work(dm);

    // scatter_to_zero will be destroyed at the end of scope, but it
    // will only decrement the reference counter incremented by
    // PetscObjectCompose below.
    auto vs = scatter_part(natural_work, start, length, 0);

    // this increments the reference counter of scatter_to_zero
    ierr = PetscObjectCompose((PetscObject)dm, "scatter_to_zero",
                              (PetscObject)(vs.scatter));
    PISM_CHK(ierr, "PetscObjectCompose");

    // this increments the reference counter of v_proc0
    ierr = PetscObjectCompose((PetscObject)dm, "v_proc0",
                              (PetscObject)vs.v);
    PISM_CHK(ierr, "PetscObjectCompose");

    VecScatterDestroy(&vs.scatter);

    // We DO NOT call VecDestroy(v_proc0): the petsc::Vec wrapper will
    // take care of this.
    return vs.v;
  }
  return v_proc0;
}
#endif
