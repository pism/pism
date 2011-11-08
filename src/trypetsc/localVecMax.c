
static char help[] = 
"\nDemonstrates that on DA-created *local* Vecs, the functions\n"
"VecMax(), VecMin(), VecNorm() do not work well because these\n"
"local Vecs are actually of type VECSEQ.  A solution used in PISM\n"
"is to add a PISMGlobalMax(), etc. after VecMax().\n\n"
"Note that the synopsis of DACreateLocalVector() is\n"
"'Creates a Seq PETSc vector that may be used with the DAXXX routines.'\n\n";

/*

compare single processor behavior:

$ ./localVecMax
global Vec has VecType = mpi
local Vec has VecType = seq
for global Vec:  max = 1.00, min = 0.00, infnorm = 1.00, twonorm = 6.51
for local  Vec:  max = 1.00, min = 0.00, infnorm = 1.00, twonorm = 6.51

to four processor behavior:

$ mpiexec -n 4 ./localVecMax
global Vec has VecType = mpi
local Vec has VecType = seq
for global Vec:  max = 1.00, min = 0.00, infnorm = 1.00, twonorm = 6.51
for local  Vec:  max = 0.50, min = 0.00, infnorm = 0.50, twonorm = 1.82

*/

#include "petscdmda.h"
#include <math.h>

extern PetscErrorCode fillWithData(DA,Vec*,Vec*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  PetscReal      maxglobal,minglobal,maxlocal,minlocal,
                 infnormglobal, twonormglobal, infnormlocal, twonormlocal;
  Vec            xglobal,xlocal;
  DA             da;

  PetscInitialize(&argc,&argv,(char *)0,help);

  /* here we follow ex5.c in creation of DA, but the same behavior
     is also seen with DA_PERIODICXY and with the tranpose used in PISM */
  ierr = DACreate2d(PETSC_COMM_WORLD, DA_NONPERIODIC,DA_STENCIL_STAR,
                    -11,-11,PETSC_DECIDE,PETSC_DECIDE, 
                    1,1,PETSC_NULL,PETSC_NULL,&da); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(da,&xglobal);CHKERRQ(ierr);
  ierr = DACreateLocalVector(da,&xlocal);CHKERRQ(ierr);

  ierr = fillWithData(da,&xglobal,&xlocal); CHKERRQ(ierr);

  // I was unable to make VecGetType() work, but direct access to Vec as a struct* works:
  ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "global Vec has VecType = %s\nlocal Vec has VecType = %s\n",
            ((PetscObject)xglobal)->type_name,((PetscObject)xlocal)->type_name);
            CHKERRQ(ierr);
  
  ierr = VecMax(xglobal,PETSC_NULL,&maxglobal); CHKERRQ(ierr);
  ierr = VecMin(xglobal,PETSC_NULL,&minglobal); CHKERRQ(ierr);
  ierr = VecMax(xlocal,PETSC_NULL,&maxlocal); CHKERRQ(ierr);
  ierr = VecMin(xlocal,PETSC_NULL,&minlocal); CHKERRQ(ierr);
  ierr = VecNorm(xglobal,NORM_INFINITY,&infnormglobal); CHKERRQ(ierr);
  ierr = VecNorm(xlocal,NORM_INFINITY,&infnormlocal); CHKERRQ(ierr);
  ierr = VecNorm(xglobal,NORM_2,&twonormglobal); CHKERRQ(ierr);
  ierr = VecNorm(xlocal,NORM_2,&twonormlocal); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "for global Vec:  max = %.2f, min = %.2f, infnorm = %.2f, twonorm = %.2f\n",
            maxglobal, minglobal, infnormglobal, twonormglobal);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, 
            "for local  Vec:  max = %.2f, min = %.2f, infnorm = %.2f, twonorm = %.2f\n",
            maxlocal, minlocal, infnormlocal, twonormlocal);CHKERRQ(ierr);

  ierr = VecDestroy(xglobal);CHKERRQ(ierr);
  ierr = VecDestroy(xlocal);CHKERRQ(ierr);      

  ierr = PetscFinalize();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode fillWithData(DA da, Vec *f, Vec *g) {
  PetscErrorCode ierr;
  PetscScalar **ff, **gg;
  PetscScalar dx,dy,xx,yy;
  PetscInt Mx,My,xs,ys,xm,ym,i,j;
  
  ierr = DAGetInfo(da,PETSC_IGNORE,&Mx,&My,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
            PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&xs,&ys,PETSC_NULL,&xm,&ym,PETSC_NULL);CHKERRQ(ierr);
  dx     = 1.0/(PetscReal)(Mx-1);
  dy     = 1.0/(PetscReal)(My-1);

  ierr = DAVecGetArray(da, *f, &ff); CHKERRQ(ierr);
  ierr = DAVecGetArray(da, *g, &gg); CHKERRQ(ierr);

  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      xx = (PetscReal)(i) * dx;
      yy = (PetscReal)(j) * dy;
      ff[j][i] = xx;
      gg[j][i] = xx;
    }
  }
  
  ierr = DAVecRestoreArray(da, *f, &ff); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da, *g, &gg); CHKERRQ(ierr);

  return 0; 
} 

