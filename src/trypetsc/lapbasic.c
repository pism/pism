
static char help[] = 
   "Solves a 2-dimensional periodic boundary condition, Poisson-like equation\n"
   "using only Vec, Mat, DA, and KSP concepts.  (No DMMG or SNES.)\n"
   "Good basic illustration of use of DAGetMatrix and MatValuesSetStencil.\n"
   "Relates Vec and Mat ownership ranges.  Based very loosely on KSP examples\n"
   "ex23.c and ex32.c.\n\n";

/*!
The problem being solved here is 
  u_xx + u_yy + u = f(x,y)    on  [0,1] x [0,1]
with periodic boundary conditions.

We set
  f(x,y) = (1 - 20 pi^2) sin(2 pi x) cos(4 pi y)
so that
  u(x,y) = sin(2 pi x) cos(4 pi y)
is the exact solution.
 */
 
/* runs:

10^7 var matrix problem done pretty quickly:
   mpiexec -n 2 ./lapbasic -m 3000 -pc_type jacobi

compare graphical views of assembly:
   mpiexec -n 4 ./lapbasic -mat_view_draw -draw_pause 5 -display :0
   mpiexec -n 4 ./lapbasic -pc_type jacobi -mat_view_draw -draw_pause 5 -display :0
   mpiexec -n 4 ./poisson -mat_view_draw -draw_pause 5 -display :0

*/

#include "petscksp.h"
#include "petscda.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscInitialize(&argc,&args,(char *)0,help);

  PetscInt       m = 10;  /* default number of rows and columns IN GRID,
                             but matrix is (m^2) x (m^2)                 */
  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);

  MPI_Comm    com = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* create m x m two-dimensional grid for periodic boundary condition problem */
  DA da2;
  PetscInt dof=1, stencilwidth=1;
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_STAR,
                    m,m,PETSC_DECIDE,PETSC_DECIDE, 
                    dof,stencilwidth,PETSC_NULL,PETSC_NULL,&da2); CHKERRQ(ierr);

  /* get da2-managed Vecs */
  Vec x,b,u;
  ierr = DACreateGlobalVector(da2,&x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&u); CHKERRQ(ierr);

  Mat A;
  ierr = DAGetMatrix(da2, MATMPIAIJ, &A); CHKERRQ(ierr);

  /* alternative call below is not quite same as result from DAGetMatrix(),
     because of nonzero allocation; the Mat ownership ranges are same       */
  /* ierr = MatCreateMPIAIJ(com, mlocal, mlocal, m*m, m*m,
                            5, PETSC_NULL, 4, PETSC_NULL,
                            &A); CHKERRQ(ierr)              */

  ierr = MatSetFromOptions(A);CHKERRQ(ierr);

  /* report on ownership range */
  PetscInt rstart,rend,mlocal;
  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&mlocal);CHKERRQ(ierr);
  PetscInt A_rstart,A_rend;
  ierr = MatGetOwnershipRange(A,&A_rstart,&A_rend);CHKERRQ(ierr);
  if ((rstart != A_rstart) || (rend != A_rend)) {
    ierr = PetscPrintf(com,
      "Vec and Mat ownership ranges different!!!  ending ...\n"); CHKERRQ(ierr);
    PetscEnd();
  } else {
    ierr = PetscSynchronizedPrintf(com,
      "rank=%d has Vec and Mat ownership:   mlocal=%d, rstart=%d, rend=%d\n",
      rank,mlocal,rstart,rend); CHKERRQ(ierr);
  }
  PetscSynchronizedFlush(com);

  /* get local part of grid */
  PetscInt  xm,ym,xs,ys;
  DAGetCorners(da2,&xs,&ys,0,&xm,&ym,0);

  /* report on local part of grid */
  ierr = PetscSynchronizedPrintf(com,
    "rank=%d has da2-managed-Vec local ranges:   xs=%d, xm=%d, ys=%d, ym=%d\n",
    rank,xs,xm,ys,ym); CHKERRQ(ierr);
  PetscSynchronizedFlush(com);

  /* set up linear system */
  PetscScalar **barr, **uarr;  /* RHS and exact soln, resp. */
  DAVecGetArray(da2, b, &barr);
  DAVecGetArray(da2, u, &uarr);

  PetscScalar dx = 1.0/(double)m, dy = dx, pi = 3.14159265358979;
  PetscScalar xi,yj;

  PetscInt    diag=0,north=1,east=2,south=3,west=4;
  PetscScalar vals[5] = {-4.0 + dx * dx, 1.0, 1.0, 1.0, 1.0};
  MatStencil  row, col[5];  /* these are not "stencils" at all, but local grid
                               to global indices helpers */

  PetscInt  i,j,num;
  for (j=ys; j<ys+ym; j++)  {
    for(i=xs; i<xs+xm; i++) {
    
      /* entries of matrix A */
      row.i = i; row.j = j; row.c = 0;  /* dof = 1 so first component; 
                                           note row.k is for 3d DAs    */
      for (num=0; num<5; num++)   col[num].c = 0;
      /* set diag first, then go through stencil neighbors */
      col[diag].i  = i;   col[diag].j  = j;
      col[north].i = i;   col[north].j = j+1; 
      col[east].i  = i+1; col[east].j  = j; 
      col[south].i = i;   col[south].j = j-1; 
      col[west].i  = i-1; col[west].j  = j; 
      ierr = MatSetValuesStencil(A,1,&row,5,col,vals,INSERT_VALUES); CHKERRQ(ierr);

      /* entries of vectors: exact solution u and right-hand-side b */
      xi = (double)i * dx;
      yj = (double)j * dy;
      uarr[j][i] = sin(2.0 * pi * xi) * cos(4.0 * pi * yj); 
      barr[j][i] = (1.0 - 20.0 * pi * pi) * uarr[j][i];
      barr[j][i] *= dx * dx;

    }
  }

  DAVecRestoreArray(da2, b, &barr);
  DAVecRestoreArray(da2, u, &uarr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u); CHKERRQ(ierr);

  /* uncomment for dense view; -mat_view default format is good enough for most purposes
  PetscViewer viewer;
  PetscViewerCreate(com, &viewer);
  PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE);
  MatView(A,viewer);
  PetscViewerDestroy(viewer);
  */

  /* setup solver context now that Mat and Vec are assembled */
  KSP ksp;
  ierr = KSPCreate(com,&ksp);CHKERRQ(ierr);

  /* Set "operators".  Here the matrix that defines the linear system
     also serves as the preconditioning matrix.  But we do not assert a
     relationship between their nonzero patterns.(???)                     */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /* Following is optional; parameters could be set at runtime.  */
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    CHKERRQ(ierr);

  /*  Set runtime options, e.g.,
    -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
 
  /*  Solve linear system  */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr); 

  /* Compute and report the error (and the iteration count and reason). */
  PetscScalar norminf, normtwo, neg_one=-1.0;
  PetscInt    its;
  KSPConvergedReason reason;
  ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);              // x = x - u
  ierr = VecNorm(x,NORM_INFINITY,&norminf);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&normtwo);CHKERRQ(ierr);           // discrete norm
  normtwo *= dx * dy;                                         // integral norm
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp,&reason); CHKERRQ(ierr);
  ierr = PetscPrintf(com,
     "Error norms  ||err||_inf = %.3e, ||err||_2 = %.3e;\n"
     "Iterations = %d;  Reason = %d\n",
     norminf, normtwo, its, (int) reason);CHKERRQ(ierr);

  /* destroy */
  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(u);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = DADestroy(da2);CHKERRQ(ierr);

  /* Always call PetscFinalize() before exiting a program. */
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

