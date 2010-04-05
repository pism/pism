
static char help[] = 
   "Solves a 2-spatial-dimension AND 2-degree-of-freedom (vectorized),\n"
   "periodic boundary condition problem.  A vector-Poisson-like equation.\n"
   "Uses only Vec, Mat, DA, and KSP concepts.  (No DMMG or SNES.)\n"
   "Makes good basic illustration of use of MatValuesSetStencil in dof>1 case.\n"
   "Shows Vec and Mat ownership ranges.  Compare lapbasic.c\n\n";

/*!
The problem being solved here is 
  u_xx + u_yy +   u - 2 v = f(x,y)   
  v_xx + v_yy - 2 u +   v = g(x,y)   
On (x,y) \in [0,1] x [0,1] with periodic boundary conditions.

We set (manufacture) an exact solution:
  u(x,y) = cos(4 pi y),
  v(x,y) = sin(2 pi x),
so
  f(x,y) = (1 - 16 pi^2) u(x,y) - 2 v(x,y)
  g(x,y) = (1 - 4 pi^2) v(x,y) - 2 u(x,y)

Runs:

compare graphical views of assembly:
   mpiexec -n 4 ./lapdof2 -mat_view_draw -draw_pause 5 -display :0

Discretization:  We require dx=dy.  The first equation becomes

u_i+1,j + u_i-1,j + u_i,j+1 + u_i,j-1 + (-4+dx^2) u_ij - 2 dx^2 v_ij
    = dx^2 f_ij

The second equation is similar.  We use a dof=2 DA for the single Vec x which
stores both u and v.  The exact solution Vec uv and the right-hand-side Vec b
have the same grid meaning.  Noting that "m" is an input and runtime option, the 
finite difference grid is m by m so x,uv,b are in a 2m^2 dimension vector space.
The storage order is something like
  x_0    = u_00
  x_1    = v_00
  x_2    = u_01
  x_3    = v_01
  x_4    = u_02
  ...
  x_2m   = u_10
  x_2m+1 = v_10
  x_2m+2 = u_11
  ...
  x_2m^2-1 = v_m-1,m-1
but basically we don't need to know that.  When we insert into u,b we 

The same DA which manages x also gives us the Mat A, with its automatically-
determined nonzero pattern.  The row and column indices of A have the same
i,j-grid-to-global indexing issues which apply to x,uv,b.  The "MatStencil"
structs are used to insert values into A, and this example shows their use when
dof=2.

*An important point is to not think of these as *stencils* in the numerical
analysis sense, but merely as little tools helping the translation of grid
coordinates i,j to the global row and column in A.*
*/

#include "petscksp.h"
#include "petscda.h"

/* the dof=2 Vecs are m-length arrays of these pairs */
typedef struct {
  PetscScalar u, v;
} pair;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscInitialize(&argc,&args,(char *)0,help);

  PetscInt       m = 10;  /* default number of rows and columns IN THE GRID,
                             while matrix is (2 m^2) x (2 m^2)                 */
  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);

  MPI_Comm    com = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* create m x m two-dimensional grid for periodic boundary condition problem */
  DA da22;

  PetscInt dof=2, stencilwidth=1;
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_STAR,
                    m,m,PETSC_DECIDE,PETSC_DECIDE, 
                    dof,stencilwidth,PETSC_NULL,PETSC_NULL,&da22); CHKERRQ(ierr);

  /* get da22-managed Vecs */
  Vec x,b,uv;
  ierr = DACreateGlobalVector(da22,&x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&uv); CHKERRQ(ierr);

  Mat A;
  ierr = DAGetMatrix(da22, MATMPIAIJ, &A); CHKERRQ(ierr);
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
  DAGetCorners(da22,&xs,&ys,0,&xm,&ym,0);

  /* report on local part of grid */
  ierr = PetscSynchronizedPrintf(com,
    "rank=%d has da22-managed-Vec local ranges:   xs=%d, xm=%d, ys=%d, ym=%d\n",
    rank,xs,xm,ys,ym); CHKERRQ(ierr);
  PetscSynchronizedFlush(com);

  /* set up linear system */
  pair **barr, **uvarr;  /* RHS and exact soln, resp. */
  DAVecGetArray(da22, b, &barr);
  DAVecGetArray(da22, uv, &uvarr);

  PetscScalar dx = 1.0/(double)m, dy = dx, pi = 3.14159265358979;
  PetscScalar xi,yj;

  MatStencil  row, col[6];  /* these are not "stencils" at all, but local grid
                                   to global indices helpers */
  PetscScalar dxsqr = dx * dx;
  PetscInt diag=0,    nort=1, east=2, sout=3, west=4, other=5;
  PetscScalar vals[10] = {-4.0+dxsqr, 1.0, 1.0, 1.0, 1.0, -2.0*dxsqr};
  
  PetscInt  i,j;
  for (j=ys; j<ys+ym; j++)  {
    for(i=xs; i<xs+xm; i++) {
    
      /* illustrations of dof=2 use of MatStencil! */
      
      /* pattern in 1st component equation:    u_xx + u_yy + u - 2 v = f   */
      row.i = i; row.j = j; row.c = 0;
      col[diag].i = i;   col[diag].j = j;    col[diag].c = 0;
      col[nort].i = i;   col[nort].j = j+1;  col[nort].c = 0;
      col[east].i = i+1; col[east].j = j;    col[east].c = 0;
      col[sout].i = i;   col[sout].j = j-1;  col[sout].c = 0; 
      col[west].i = i-1; col[west].j = j;    col[west].c = 0;
      col[other].i = i;  col[other].j = j;   col[other].c = 1;

      /* set entries of matrix A for these two equations */
      ierr = MatSetValuesStencil(A,1,&row,6,col,vals,INSERT_VALUES); CHKERRQ(ierr);

      /* pattern in 2nd component equation:   v_xx + v_yy - 2 u + v = g   */
      row.i = i; row.j = j; row.c = 1;
      col[diag].i = i;   col[diag].j = j;    col[diag].c = 1;
      col[nort].i = i;   col[nort].j = j+1;  col[nort].c = 1;
      col[east].i = i+1; col[east].j = j;    col[east].c = 1;
      col[sout].i = i;   col[sout].j = j-1;  col[sout].c = 1; 
      col[west].i = i-1; col[west].j = j;    col[west].c = 1;
      col[other].i = i;  col[other].j = j;   col[other].c = 0;

      /* set entries of matrix A for these two equations */
      ierr = MatSetValuesStencil(A,1,&row,6,col,vals,INSERT_VALUES); CHKERRQ(ierr);

      /* entries of vectors: exact solution u and right-hand-side b */
      xi = (double)i * dx;
      yj = (double)j * dy;
      uvarr[j][i].u = cos(4.0 * pi * yj); 
      uvarr[j][i].v = sin(2.0 * pi * xi); 
      barr[j][i].u  = (1.0 - 16.0 * pi * pi) * uvarr[j][i].u - 2.0 * uvarr[j][i].v;
      barr[j][i].u *= dx * dx;
      barr[j][i].v  = (1.0 - 4.0 * pi * pi) * uvarr[j][i].v - 2.0 * uvarr[j][i].u;
      barr[j][i].v *= dx * dx;

    }
  }

  DAVecRestoreArray(da22, b, &barr);
  DAVecRestoreArray(da22, uv, &uvarr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(uv); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(uv); CHKERRQ(ierr);

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
  ierr = VecAXPY(x,neg_one,uv);CHKERRQ(ierr);              // x = x - u
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
  ierr = VecDestroy(uv);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = DADestroy(da22);CHKERRQ(ierr);

  /* Always call PetscFinalize() before exiting a program. */
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

