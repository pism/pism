
static char help[] = "\nSolves Poisson-like problem.\n\n"
"Greatly simplified modification of code src/snes/examples/tutorials/ex5.c\n"
"in PETSc src tree, which solves Bratu nonlinear problem.\n\n"
"Here we treat a linear problem like a generic nonlinear problem.\n\n";

/* ------------------------------------------------------------------------

    Solves
  
            - epsilon Laplacian u + f(x,y) u = g(x,y),  0 < x,y < 1,
  
    with boundary conditions
   
             u = 0  for  x = 0, x = 1, y = 0, y = 1.
  
    A finite difference approximation with the usual 5-point stencil
    is used to discretize the boundary value problem to obtain a nonlinear 
    system of equations.

Program usage:  mpiexec -np <procs> poisson [-help] [all PETSc options]
e.g.:

   ./poisson -draw_pause 5 -da_grid_x 20 -da_grid_y 30 -snes_rtol 1e-2

   mpiexec -np 4 ./poisson -log_summary

   mpiexec -np 2 ./poisson -draw_pause 5 -da_grid_x 200 -da_grid_y 200 -snes_rtol 1e-20 -display :0

  ------------------------------------------------------------------------- */

/* DEMONSTRATION OF CONVERGENCE:

 $ ./poisson -da_grid_x 10 -da_grid_y 10 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 5.811e-16
Max value of abs(u_NUM-u_EXACT) = 4.522e-02
 $ ./poisson -da_grid_x 20 -da_grid_y 20 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 3.763e-16
Max value of abs(u_NUM-u_EXACT) = 1.095e-02
 $ ./poisson -da_grid_x 40 -da_grid_y 40 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 5.828e-16
Max value of abs(u_NUM-u_EXACT) = 2.590e-03
 $ ./poisson -da_grid_x 80 -da_grid_y 80 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 4.247e-16
Max value of abs(u_NUM-u_EXACT) = 6.304e-04
 $ ./poisson -da_grid_x 160 -da_grid_y 160 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 3.951e-16
Max value of abs(u_NUM-u_EXACT) = 1.558e-04
 $ ./poisson -da_grid_x 320 -da_grid_y 320 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 4.151e-16
Max value of abs(u_NUM-u_EXACT) = 3.870e-05
 $ mpiexec -np 2 ./poisson -da_grid_x 640 -da_grid_y 640 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 4.206e-16
Max value of abs(u_NUM-u_EXACT) = 9.644e-06

RELATED DEMOS OF GOOD BEHAVIOR:
 $ mpiexec -np 2 ./poisson -da_grid_x 320 -da_grid_y 320 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 4.077e-16
Max value of abs(u_NUM-u_EXACT) = 3.870e-05
 $ mpiexec -np 2 ./poisson -da_grid_x 320 -da_grid_y 320 -snes_rtol 1.0e-2
Number of Newton iterations = 1
Max value of abs(residual) = 4.561e-07
Max value of abs(u_NUM-u_EXACT) = 3.935e-05
 $ mpiexec -np 2 ./poisson -da_grid_x 320 -da_grid_y 207 -snes_rtol 1.0e-20
Number of Newton iterations = 3
Max value of abs(residual) = 4.477e-16
Max value of abs(u_NUM-u_EXACT) = 9.213e-05

 */

#include "petscda.h"
#include "petscsnes.h"
#include <math.h>

/* 
   User-defined application context - contains data needed by the 
   application-provided call-back routines, FormJacobianLocal() and
   FormFunctionLocal().
*/
typedef struct {
   DA          da;             /* distributed array data structure */
   PetscReal   epsilon;
   Vec         f;              /* = f(x,y) in above */
   Vec         g;              /* = g(x,y) in above */
   Vec         uexact;         /* = u(x,y), exact value, in above */
} AppCtx;

/*    User-defined routines  */
extern PetscErrorCode fillPoissonData(AppCtx*);
extern PetscErrorCode FormFunctionLocal(DALocalInfo*,PetscScalar**,PetscScalar**,AppCtx*);
extern PetscErrorCode FormJacobianLocal(DALocalInfo*,PetscScalar**,Mat,AppCtx*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    x,r;                  /* solution, residual vectors */
  Mat                    A,J;                    /* Jacobian matrix */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               its;                  /* iterations for convergence */

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,-4,-4,PETSC_DECIDE,PETSC_DECIDE,
                    1,1,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.f);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.g);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.uexact);CHKERRQ(ierr);

  /* main added content: */
  ierr = fillPoissonData(&user); CHKERRQ(ierr);
  
  ierr = DAGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
  
  A    = J;

  ierr = SNESSetJacobian(snes,A,J,SNESDAComputeJacobian,&user);CHKERRQ(ierr); // default

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)FormFunctionLocal);CHKERRQ(ierr);
  ierr = DASetLocalJacobian(user.da,(DALocalFunction1)FormJacobianLocal);CHKERRQ(ierr); 

  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = VecSet(x,0.0); CHKERRQ(ierr); /* REPLACES in ex5.c: FormInitialGuess(&user,x); */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSolve(snes,PETSC_NULL,x);CHKERRQ(ierr); 

  /* feedback */
  PetscReal resnorm, errnorm;
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of Newton iterations = %D\n",its);
            CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_INFINITY,&resnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Max value of abs(residual) = %9.3e\n",resnorm);
            CHKERRQ(ierr);

/* OPTIONAL viewer stuff:
  PetscViewer viewer;
  ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,"solution x",
            PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,&viewer); CHKERRQ(ierr);
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
*/
  ierr = VecAXPY(x,-1.0,user.uexact); CHKERRQ(ierr);
/* final OPTIONAL viewer stuff:
  ierr = VecView(x,viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
*/

  ierr = VecNorm(x,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Max value of abs(u_NUM-u_EXACT) = %9.3e\n",errnorm);CHKERRQ(ierr);
  
  if (A != J) {
    ierr = MatDestroy(A);CHKERRQ(ierr);
  }
  ierr = MatDestroy(J);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);      
  ierr = VecDestroy(user.f);CHKERRQ(ierr);      
  ierr = VecDestroy(user.g);CHKERRQ(ierr);      
  ierr = VecDestroy(user.uexact);CHKERRQ(ierr);      
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
                   PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
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

      ff[j][i] = 1000.0 * xx * yy;
      uuex[j][i] = sin(pi * xx) * sin(3.0 * pi * yy);
      gg[j][i] = (user->epsilon * 10.0 * pi * pi + ff[j][i]) * uuex[j][i];

/* EASIER CASE:
      ff[j][i] = 0.0;
      uuex[j][i] = sin(pi * xx) * sin(pi * yy);
      gg[j][i] = (user->epsilon * 2.0 * pi * pi) * uuex[j][i];
*/
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
        v[0] = -scyy;                                           col[0].j = j - 1; col[0].i = i;
        v[1] = -scxx;                                           col[1].j = j;     col[1].i = i-1;
        v[2] = eps * 2.0 * (scxx + scyy) + sc * ff[j][i];       col[2].j = row.j; col[2].i = row.i;
        v[3] = -scxx;                                           col[3].j = j;     col[3].i = i+1;
        v[4] = -scyy;                                           col[4].j = j + 1; col[4].i = i;
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

