
static char help[] = "\nSolves Poisson-like problem.\n\n"
"Greatly simplified modification of code src/snes/examples/tutorials/ex5.c\n"
"in PETSc src tree, which solves Bratu nonlinear problem.\n\n"
"Here we treat a linear problem like a generic nonlinear problem.\n\n";

/* Solves

            - epsilon Laplacian u + f(x,y) u = g(x,y),  0 < x,y < 1,
  
with boundary conditions
   
             u = 0  for  x = 0, x = 1, y = 0, y = 1.
  
A finite difference approximation with the usual 5-point stencil is used 
to discretize the boundary value problem to obtain a nonlinear system of
equations.

Program usage:  mpiexec -np <procs> poisson [-help] [all PETSc options]
e.g.:

   ./poisson -da_grid_x 20 -da_grid_y 30 -snes_rtol 1e-2

   mpiexec -np 2 ./poisson -dodraw -draw_pause 5 -display :0

--------------------------------------------
DEMONSTRATION OF CONVERGENCE:

$ for M in 20 40 80 160 320 640; do 
>   mpiexec -np 2 ./poisson -da_grid_x $M -da_grid_y $M -snes_rtol 1e-15; done
 20 x  20 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 3.970e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 1.095e-02
 40 x  40 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 4.110e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 2.590e-03
 80 x  80 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 5.942e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 6.304e-04
160 x 160 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 3.856e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 1.558e-04
320 x 320 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 4.125e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 3.870e-05
640 x 640 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 4.206e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 9.644e-06

RELATED DEMOS OF GOOD BEHAVIOR (PARAMETER CHANGES):
$ ./poisson -da_grid_x 320 -da_grid_y 207 -snes_rtol 1.0e-15  # NONEQUAL GRID
320 x 207 grid: number of Newton iterations = 3
                max abs(residual) = |residual|_infty = 4.477e-16
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 9.213e-05
$ mpiexec -np 2 ./poisson -da_grid_x 320 -da_grid_y 320 -snes_rtol 1.0e-2 # WEAK RTOL
320 x 320 grid: number of Newton iterations = 1
                max abs(residual) = |residual|_infty = 4.561e-07
                max abs(u_NUM-u_EXACT) = |u_NUM - u_EXACT|_infty = 3.935e-05

 */

#include "petscda.h"
#include "petscsnes.h"
#include <math.h>

/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal(), and by fillPoissonData().
*/
typedef struct {
   DA          da;             /* must be first in struct */
   PetscReal   epsilon;
   Vec         f;              /* = f(x,y) in PDE */
   Vec         g;              /* = g(x,y) in PDE */
   Vec         uexact;         /* = u(x,y), exact value, in PDE */
} AppCtx;

extern PetscErrorCode fillPoissonData(AppCtx*);
extern PetscErrorCode FormFunctionLocal(DALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
extern PetscErrorCode FormJacobianLocal(DALocalInfo*,PetscScalar**,Mat,AppCtx*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    x,r;                  /* solution, residual vectors */
  Mat                    J;                    /* Jacobian matrix */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               its;                  /* iterations for convergence */

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* in ex5.c it has "DA_NONPERIODIC", but as in IceGrid works fine too: */
  ierr = DACreate2d(PETSC_COMM_WORLD, DA_XYPERIODIC,DA_STENCIL_STAR,
                    -10,-10,PETSC_DECIDE,PETSC_DECIDE, 
                    1,1,PETSC_NULL,PETSC_NULL,&user.da); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.f);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.g);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.uexact);CHKERRQ(ierr);

  /* main added content re ex5: nontrivial coeffs */
  ierr = fillPoissonData(&user); CHKERRQ(ierr);
  
  ierr = DAGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
  
  /* use default method of Jacobian eval (i.e. uses FormJacobianLocal because of
     DASetLocalJacobian() below); also preconditioner is same as Jacobian;
     compare different approaches here in ex5.c */
  ierr = SNESSetJacobian(snes,J,J,SNESDAComputeJacobian,&user);CHKERRQ(ierr); // default

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)FormFunctionLocal);CHKERRQ(ierr);

  ierr = DASetLocalJacobian(user.da,(DALocalFunction1)FormJacobianLocal);CHKERRQ(ierr); 

  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = VecSet(x,0.0); CHKERRQ(ierr); /* compare in ex5.c: FormInitialGuess(&user,x); */

  /* solve linear system (by nonlinear means)! */
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr); 

  /* feedback: grid, numerical method convergence info, and error */
  PetscInt Mx,My;
  PetscReal resnorm, uerrnorm;
  Vec uerr;
  ierr = DAGetInfo(user.da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr); 
  ierr = VecNorm(r,NORM_INFINITY,&resnorm); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&uerr);CHKERRQ(ierr);
  ierr = VecWAXPY(uerr,-1.0,x,user.uexact); CHKERRQ(ierr); /* uerr = -x + user.uexact */
  ierr = VecNorm(uerr,NORM_INFINITY,&uerrnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
            "%3d x %3d grid: number of Newton iterations = %d\n",Mx,My,its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
            "                max abs(residual) = |residual|_infty = %9.3e\n",resnorm);
            CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
            "                max abs(u_NUM-u_EXACT) = |u_NUM-u_EXACT|_infty = %9.3e\n",
            uerrnorm);CHKERRQ(ierr);

  /* optionally, draw viewers; all fields appear in same viewer here;
     use -draw_pause N to see for N secs each, e.g. -dodraw -draw_pause 3 */
  PetscTruth dodraw = PETSC_FALSE;
  ierr = PetscOptionsGetTruth(PETSC_NULL,"-dodraw",&dodraw,0); CHKERRQ(ierr);
  if (dodraw) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
            "  -dodraw seen; showing six fields in viewer ...\n");CHKERRQ(ierr);
    PetscViewer viewer;
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,"solution x",
              PETSC_DECIDE,PETSC_DECIDE,500,500,&viewer); CHKERRQ(ierr);
    /* bug in 2.3.3: this is mechanism for titles which works ... */
    PetscDraw draw;
    ierr = PetscViewerDrawGetDraw(viewer,0,&draw); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"num solution u(x,y)"); CHKERRQ(ierr);
    ierr = VecView(x,viewer); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"residual r(x,y)"); CHKERRQ(ierr);
    ierr = VecView(r,viewer); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"f(x,y)"); CHKERRQ(ierr);
    ierr = VecView(user.f,viewer); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"g(x,y)"); CHKERRQ(ierr);
    ierr = VecView(user.g,viewer); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"exact solution u(x,y)"); CHKERRQ(ierr);
    ierr = VecView(user.uexact,viewer); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw,"error (u_NUM-u_EXACT)"); CHKERRQ(ierr);
    ierr = VecView(uerr,viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
  }

  /* de-allocate */
  ierr = MatDestroy(J);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);      
  ierr = VecDestroy(user.f);CHKERRQ(ierr);      
  ierr = VecDestroy(user.g);CHKERRQ(ierr);      
  ierr = VecDestroy(user.uexact);CHKERRQ(ierr);      
  ierr = VecDestroy(uerr); CHKERRQ(ierr);
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = DADestroy(user.da);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode fillPoissonData(AppCtx* user) {
  PetscErrorCode ierr;
  PetscScalar **ff, **gg, **uuex;
  PetscScalar dx,dy,xx,yy,pi;
  PetscInt Mx,My,xs,ys,xm,ym,i,j;

  user->epsilon = 1.0;
  
  pi = 3.14159265358979;
  ierr = DAGetInfo(user->da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  dx     = 1.0/(PetscReal)(Mx-1);
  dy     = 1.0/(PetscReal)(My-1);
  ierr = DAGetCorners(user->da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);

  ierr = DAVecGetArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da, user->g, &gg); CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da, user->uexact, &uuex); CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      xx = (PetscReal)(i) * dx;
      yy = (PetscReal)(j) * dy;
      
      /* here we manufacture a solution; i.e. we choose g(x,y)
         so that the given f(x,y) and u(x,y)=u_EXACT(x,y) are
         a solution.
         Note user->uexact is only evaluated at the end, for error,
         not during numerical run. */
      ff[j][i] = 1000.0 * xx * yy;
      uuex[j][i] = sin(pi * xx) * sin(3.0 * pi * yy);
      gg[j][i] = (user->epsilon * 10.0 * pi * pi + ff[j][i]) * uuex[j][i];
    }
  }
  ierr = DAVecRestoreArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da, user->g, &gg); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da, user->uexact, &uuex); CHKERRQ(ierr);

  return 0; 
} 


/*   FormFunctionLocal - Evaluates nonlinear function, F(x). */
PetscErrorCode FormFunctionLocal(DALocalInfo *info,PetscScalar **x,PetscScalar **F,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  PetscReal      dx,dy,sc,scxx,scyy,eps;
  PetscScalar    u,neguxx,neguyy;
  PetscScalar    **ff, **gg;

  PetscFunctionBegin;

  dx     = 1.0/(PetscReal)(info->mx-1);
  dy     = 1.0/(PetscReal)(info->my-1);
  sc     = dx*dy;
  scxx   = sc/(dx*dx);
  scyy   = sc/(dy*dy);
  eps    = user->epsilon;

  ierr = DAVecGetArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da, user->g, &gg); CHKERRQ(ierr);
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        F[j][i] = x[j][i];
      } else {
        u       = x[j][i];
        neguxx  = (2.0*u - x[j][i-1] - x[j][i+1])*scxx;
        neguyy  = (2.0*u - x[j-1][i] - x[j+1][i])*scyy;
        F[j][i] = eps * (neguxx + neguyy) + sc * ff[j][i] * u - sc * gg[j][i];
      }
    }
  }
  ierr = DAVecRestoreArray(user->da, user->f, &ff); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da, user->g, &gg); CHKERRQ(ierr);

  PetscFunctionReturn(0); 
} 


/*   FormJacobianLocal - Evaluates Jacobian matrix.  */
PetscErrorCode FormJacobianLocal(DALocalInfo *info,PetscScalar **x,Mat jac,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,j;
  MatStencil     col[5],row;
  PetscScalar    v[5],dx,dy,sc,scxx,scyy,eps;
  PetscScalar    **ff;

  PetscFunctionBegin;

  dx     = 1.0/(PetscReal)(info->mx-1);
  dy     = 1.0/(PetscReal)(info->my-1);
  sc     = dx*dy;
  scxx = sc/(dx * dx);
  scyy = sc/(dy * dy);
  eps  = user->epsilon;

  ierr = DAVecGetArray(user->da, user->f, &ff); CHKERRQ(ierr);
  for (j=info->ys; j<info->ys+info->ym; j++) {
    for (i=info->xs; i<info->xs+info->xm; i++) {
      row.j = j; row.i = i;
      /* boundary points */
      if (i == 0 || j == 0 || i == info->mx-1 || j == info->my-1) {
        v[0] = 1.0;
        ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
      /* interior grid points */
        v[0] = -scyy;                                     col[0].j = j - 1; col[0].i = i;
        v[1] = -scxx;                                     col[1].j = j;     col[1].i = i-1;
        v[2] = eps * 2.0 * (scxx + scyy) + sc * ff[j][i]; col[2].j = row.j; col[2].i = row.i;
        v[3] = -scxx;                                     col[3].j = j;     col[3].i = i+1;
        v[4] = -scyy;                                     col[4].j = j + 1; col[4].i = i;
        ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(user->da, user->f, &ff); CHKERRQ(ierr);

  /* Assemble matrix, using the 2-step process */
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  /* Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error. */
  ierr = MatSetOption(jac,MAT_NEW_NONZERO_LOCATION_ERR);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

