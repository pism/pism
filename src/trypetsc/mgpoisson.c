static char help[] = "Simple Poisson equation in 2D, with Dirichlet boundary conditions:\n"
"   - Laplacian u = 1  on   0 < x,y < 1,\n"
"               u = 1  for  x = 0 | x = 1 | y = 0 | y = 1.\n"
"Uses multigrid = DMMG to solve the linear system.\n"
"\nExample based on src/ksp/ksp/examples/tutorials/ex22.c.\n"
"\nRuns:\n"
"  mpiexec -n 4 ./mgpoisson -dmmg_nlevels 5 -dmmg_view\n"
"  ./mgpoisson -dmmg_nlevels 3 -mat_view_draw -draw_pause 2\n"
"\nWeak-scaling study:\n"
"  mpiexec -n 1 ./mgpoisson -dmmg_nlevels 4 -log_summary\n"
"  mpiexec -n 4 ./mgpoisson -dmmg_nlevels 5 -log_summary\n"
"  mpiexec -n 16 ./mgpoisson -dmmg_nlevels 6 -log_summary\n"
"\n";

#include "petscksp.h"
#include "petscdmmg.h"

extern PetscErrorCode ComputeMatrix(DMMG,Mat,Mat);
extern PetscErrorCode ComputeRHS(DMMG,Vec);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  DMMG           *dmmg;
  PetscReal      norm;
  PetscInt       mx,my;
  DA             da;

  PetscInitialize(&argc,&argv,(char *)0,help);

  /* create multigrid object that defaults to 2 levels */
  ierr = DMMGCreate(PETSC_COMM_WORLD,2,PETSC_NULL,&dmmg);CHKERRQ(ierr);

  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
                    -33,-33,                      /* m,n default to 10 */
                    PETSC_DECIDE,PETSC_DECIDE,  /* let PETSc choose proc per dim */
                    1,1,PETSC_NULL,PETSC_NULL,  /* dof=1, stencil width=1 */
                    &da);CHKERRQ(ierr);  
  ierr = DMMGSetDM(dmmg,(DM)da);CHKERRQ(ierr);
  ierr = DADestroy(da);CHKERRQ(ierr);

  ierr = DMMGSetKSP(dmmg,ComputeRHS,ComputeMatrix);CHKERRQ(ierr);

  ierr = DMMGSetUp(dmmg);CHKERRQ(ierr);
  ierr = DAGetInfo(DMMGGetDA(dmmg),0,&mx,&my,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"fine grid is %d x %d points\n",mx,my);CHKERRQ(ierr);

  ierr = DMMGSolve(dmmg);CHKERRQ(ierr);

  ierr = MatMult(DMMGGetJ(dmmg),DMMGGetx(dmmg),DMMGGetr(dmmg));CHKERRQ(ierr);
  ierr = VecAXPY(DMMGGetr(dmmg),-1.0,DMMGGetRHS(dmmg));CHKERRQ(ierr);
  ierr = VecNorm(DMMGGetr(dmmg),NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"residual norm %G\n",norm);CHKERRQ(ierr);

  ierr = DMMGDestroy(dmmg);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(DMMG dmmg,Vec b)
{
  PetscErrorCode ierr;
  PetscInt       mx,my;
  PetscScalar    h;

  PetscFunctionBegin;
  ierr = DAGetInfo((DA)dmmg->dm,0,&mx,&my,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  h    = 1.0/((mx-1)*(my-1));
  ierr = VecSet(b,h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
    
#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(DMMG dmmg,Mat jac,Mat B)
{
  DA             da = (DA)dmmg->dm;
  PetscErrorCode ierr;
  PetscInt       i,j, mx,my, xm,ym, xs,ys;
  PetscScalar    hx,hy, hxdhy,hydhx, cc;  /* for scaling eqns to have O(1) coeffs */
  PetscScalar    v[5];
  MatStencil     row,col[5];

  ierr = DAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0);CHKERRQ(ierr);  
  hx = 1.0 / (PetscScalar)(mx-1);  hy = 1.0 / (PetscScalar)(my-1);
  hxdhy = hx/hy;  hydhx = hy/hx;
  cc = 2.0*(hxdhy + hydhx);
  
  ierr = DAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++){
    for(i=xs; i<xs+xm; i++){
      row.i = i; row.j = j;
      if (i==0 || j==0 || i==mx-1 || j==my-1){
        v[0] = cc;
	ierr = MatSetValuesStencil(B,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
      } else {
	v[0] = -hxdhy; col[0].i = i;   col[0].j = j-1;
	v[1] = -hydhx; col[1].i = i-1; col[1].j = j;
	v[2] = cc;     col[2].i = i;   col[2].j = j;
	v[3] = -hydhx; col[3].i = i+1; col[3].j = j;
	v[4] = -hxdhy; col[4].i = i;   col[4].j = j+1;
	ierr = MatSetValuesStencil(B,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}

