static const char help[] = 
"Solve Bodvardsson equations (Bueler interpretation) using SNES, a dof=2 Vec\n"
"holding both thickness H and velocity u, and finite differences to approximate\n"
"the mass continuity and SSA stress balance PDEs.\n\n"
"This run certainly looks like success:\n"
"  ./bodvardsson -snes_mf -da_grid_x 100 -snes_monitor -snes_monitor_solution -draw_pause 1\n"
"and for N=75,150,300,600 I get convergence in H and in u:\n"
"  ./bodvardsson -snes_fd -da_grid_x N\n"
"But note wiggles in:\n"
"  ./bodvardsson -snes_fd -da_grid_x 30 -snes_monitor -snes_monitor_solution -draw_pause 1\n"
"Pushing it:\n"
"  ./bodvardsson -snes_fd -da_grid_x 800 -pc_type lu -snes_max_funcs 100000 -snes_ls basic\n"
"See also:\n"
"  ./bodvardsson -snes_mf -da_grid_x 10 -dump\n\n"
"TO DO:  * add Picard or Jacobian\n"
"        * clean up scaling\n"
"        * reasonable initial guesses\n"
"        * try other transport schemes\n"
"        * higher order Mx-1 case for dHdx in driving stress\n";

#include "petscda.h"
#include "petscsnes.h"

#include "../exactTestN.h"

typedef struct {
  PetscScalar H, u;
} Node;

/* User-defined application context.  Used esp. by FormFunctionLocal().  */
typedef struct {
  DA          da;       /* 1d,dof=2 distributed array for soln and residual */
  DA          scalarda; /* 1d,dof=1 distributed array for parameters depending on x */
  PetscInt    Mx, xs, xm;
  PetscReal   secpera, n, rho, rhow, g;
  PetscReal   H0;      /* thickness at x=0, for Dirichlet condition on mass cont */
  PetscReal   xc;      /* location at which stress (Neumann) condition applied to SSA eqn */
  PetscReal   Txc;     /* vertically-integrated longitudinal stress at xc, for Neumann cond:
                            T = 2 H B |u_x|^{(1/n)-1} u_x  */
  PetscReal   epsilon; /* regularization of viscosity, a strain rate */
  Vec         Huexact; /* exact thickness (Huexact[i][0]) and exact velocity
                          (Huexact[i][1]) on regular grid */
  Vec         M;       /* surface mass balance on regular grid */
  Vec         beta;    /* sliding coefficient on regular grid */
  Vec         B_stag;  /* ice hardness on staggered grid*/
} AppCtx;


#undef __FUNCT__
#define __FUNCT__ "FillExactSoln"
/*  Compute the exact thickness and velocity on the regular grid. */
static PetscErrorCode FillExactSoln(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,Mx=user->Mx;
  PetscScalar    dx, x, dum1, dum2, dum3, dum4;
  Node           *Hu;

  PetscFunctionBegin;
  dx = user->xc / (PetscReal)(Mx-1);
  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DAVecGetArray(user->da,user->Huexact,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    x = dx * (PetscReal)i;  /* = x_i = distance from dome */
    ierr = exactN(x, &(Hu[i].H), &dum1, &(Hu[i].u), &dum2, &dum3, &dum4); CHKERRQ(ierr);
  }
  ierr = DAVecRestoreArray(user->da,user->Huexact,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FillDistributedParams"
/*  Compute the values of the surface mass balance, ice hardness, and sliding coeff. */
static PetscErrorCode FillDistributedParams(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       i,Mx=user->Mx;
  PetscScalar    dx, x, dum1, dum2, dum3, dum4, dum5, *M, *Bstag, *beta;

  PetscFunctionBegin;
  dx = user->xc / (PetscReal)(Mx-1);
  /* Compute regular grid exact soln and staggered-grid thickness over the
     locally-owned part of the grid */
  ierr = DAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->scalarda,user->B_stag,&Bstag);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    x = dx * (PetscReal)i;  /* = x_i = distance from dome; regular grid point */
    ierr = exactN(x, &dum1, &dum2, &dum3, &(M[i]), &dum4, &(beta[i])); CHKERRQ(ierr);
    x = x + (dx/2.0);       /* = x_{i+1/2}; staggered grid point */
    if (i < Mx-1) {
      ierr = exactN(x, &dum1, &dum2, &dum3, &dum4, &(Bstag[i]), &dum5); CHKERRQ(ierr);
    } else {
      Bstag[i] = -9999.9999;
    }
  }
  ierr = DAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->scalarda,user->B_stag,&Bstag);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ScaleVelocity"
static PetscErrorCode ScaleVelocity(AppCtx *user, Vec vHu)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Node           *Hu;

  PetscFunctionBegin;
  ierr = DAVecGetArray(user->da,vHu,&Hu);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    Hu[i].u *= user->secpera;
  }
  ierr = DAVecRestoreArray(user->da,vHu,&Hu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MaxErrorNorms"
static PetscErrorCode MaxErrorNorms(AppCtx *user, Vec firstHu, Vec secondHu,
                                    PetscReal *normH, PetscReal *normu)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscReal      lnormH, lnormu;  /* norms of local portions */
  Node           *Hu1, *Hu2;

  PetscFunctionBegin;
  lnormH = -1.0;  lnormu = -1.0;
  ierr = DAVecGetArray(user->da,firstHu,&Hu1);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,secondHu,&Hu2);CHKERRQ(ierr);
  for (i = user->xs; i < user->xs + user->xm; i++) {
    lnormH = PetscMax(lnormH, PetscAbsScalar(Hu1[i].H-Hu2[i].H));
    lnormu = PetscMax(lnormu, PetscAbsScalar(Hu1[i].u-Hu2[i].u));
  }
  ierr = DAVecRestoreArray(user->da,firstHu,&Hu1);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,secondHu,&Hu2);CHKERRQ(ierr);
  ierr = PetscGlobalMax(&lnormH,normH,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = PetscGlobalMax(&lnormu,normu,PETSC_COMM_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "GetFSR"
/* define a power of the strain rate:   F  \approx  |u_x|^q u_x   */
static inline PetscScalar GetFSR(PetscScalar dx, PetscScalar eps, PetscScalar n, 
                                 PetscScalar ul, PetscScalar ur) {
  PetscScalar dudx = (ur - ul) / dx, 
              q    = (1.0 / n) - 1.0;
  return PetscPowScalar(dudx * dudx + eps * eps, q / 2.0) * dudx;
}


#undef __FUNCT__
#define __FUNCT__ "FormFunctionLocal"
static PetscErrorCode BodFunctionLocal(DALocalInfo *info, Node *Hu, Node *f, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscReal      dx, rg, spa;
  PetscScalar    *M, *Bstag, *beta,
                 duHl, ul, u, ur, dHdx, Fl, Fr, Tl, Tr,
                 scH, scu, scuH, scstress;   /* ad hoc scaling coeffs */
  PetscInt       i, Mx = info->mx;
  Vec            locBstag;

  PetscFunctionBegin;

  /* we need stencil width on Bstag (but not for M, beta) */
  ierr = DAGetLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);  /* do NOT destroy it */
  ierr = DAGlobalToLocalBegin(user->scalarda,user->B_stag,INSERT_VALUES,locBstag); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(user->scalarda,user->B_stag,INSERT_VALUES,locBstag); CHKERRQ(ierr);

  dx = user->xc / ((PetscReal)Mx - 1.0);
  rg = user->rho * user->g;
  spa = user->secpera;

  scH = 1.0 / user->H0;
  scu = spa / 100.0;
  scuH = spa / user->H0;
  scstress = 1.0 / (rg * user->H0 * dx * 0.001);

  /* note:  we use Dirichlet conditions as values for neighbor, which tends to
            symmetrize FD- or MF-computed Jacobian */
  ierr = DAVecGetArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);
  for (i = info->xs; i < info->xs + info->xm; i++) {
  
    /* MASS CONT */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].H = Hu[0].H - user->H0;
      f[0].H *= scH;
    } else {
      duHl = (Hu[i].u/spa) * Hu[i].H - ( (i == 1) ? 0.0 : (Hu[i-1].u/spa) * Hu[i-1].H );
      f[i].H = dx * M[i] - duHl;  /* upwind; difference  (uH)  only to left */
      f[i].H *= scuH;
    }

    /* SSA */
    if (i == 0) {
      /* residual at left-most point is Dirichlet cond. */
      f[0].u = (Hu[0].u/spa) - 0.0;
      f[0].u *= scu;
    } else {
      /* residual: SSA eqn */
      /* consecutive values of u */
      ul = (i == 1) ? 0.0 : (Hu[i-1].u/spa);
      u  = (Hu[i].u/spa);
      ur = (i == Mx-1) ? -1.1e30 : (Hu[i+1].u/spa);
      /* surface slope */
      if (i == 1) { 
        dHdx  = (Hu[i+1].H - user->H0) / (2.0 * dx);
      } else if (i == Mx-1) {
        dHdx  = (Hu[i].H - Hu[i-1].H) / dx; /* IMPROVABLE */
      } else { /* generic case */
        dHdx  = (Hu[i+1].H - Hu[i-1].H) / (2.0 * dx);
      }
      /* vertically-integrated longitudinal stress */
      Fl = GetFSR(dx,user->epsilon,user->n, ul,u);
      if (i == Mx-1) {
        Tl = 2.0 * (Hu[i-1].H + Hu[i].H) * Bstag[i-1] * Fl;
        Tr = 2.0 * user->Txc;
      } else {
        Fr = GetFSR(dx,user->epsilon,user->n, u,ur);
        Tl = (Hu[i-1].H + Hu[i].H) * Bstag[i-1] * Fl;
        Tr = (Hu[i].H + Hu[i+1].H) * Bstag[i] * Fr;        
      }

      f[i].u = (Tr - Tl) - dx * beta[i] * u - dx * rg * Hu[i].H * dHdx;
      f[i].u *= scstress;

    }
  }
  ierr = DAVecRestoreArray(user->scalarda,locBstag,&Bstag);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->scalarda,user->M,&M);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->scalarda,user->beta,&beta);CHKERRQ(ierr);

  ierr = DARestoreLocalVector(user->scalarda,&locBstag);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/* TODO:
static PetscErrorCode BodJacobianMatrixLocal(DALocalInfo*,Node*,Mat,AppCtx*); */


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode         ierr;

  SNES                   snes;                 /* nonlinear solver */
  Vec                    Hu,r;                 /* solution, residual vectors */
  Mat                    J;                    /* Jacobian matrix */
  AppCtx                 user;                 /* user-defined work context */
  PetscInt               its;                  /* iteration count, num of pts */
  PetscReal              errHinf, erruinf,     /* max norm of numerical error */
                         tmp1, tmp2, tmp3, tmp4, tmp5;
  PetscTruth             eps_set = PETSC_FALSE, dump = PETSC_FALSE;
  SNESConvergedReason    reason;               /* Check convergence */

  PetscInitialize(&argc,&argv,(char *)0,help);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "BODVARDSSON solves for thickness and velocity in 1D, steady ice stream\n"
    "  [run with -help for info and options]\n");CHKERRQ(ierr);

  user.n       = 3.0;          /* Glen flow law exponent */
  user.secpera = 31556926.0;
  user.rho     = 910.0;        /* kg m^-3 */
  user.rhow    = 1028.0;       /* kg m^-3 */
  user.g       = 9.81;         /* m s^-2 */
  
  /* ask Test N for its parameters, but only those we need to solve */
  ierr = params_exactN(&(user.H0), &tmp1, &(user.xc), &tmp2, &tmp3, &tmp4, &tmp5, 
                       &(user.Txc)); CHKERRQ(ierr);

  user.epsilon = (1.0 / user.secpera) / user.xc; /* regularize using strain rate
                                                    of 1/xc = ?.?e-6 per year */

  ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  geometry: H0 = %.2f m, xc = %.2f km, Txc = %.5e Pa m\n",
      user.H0, user.xc/1000.0, user.Txc);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,
      "bodvardsson program options",__FILE__);CHKERRQ(ierr);
  {
    ierr = PetscOptionsTruth("-bod_dump",
      "dump out exact and approximate solution and residual, as asci","",
      dump,&dump,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-bod_epsilon","regularization (a strain rate in units of 1/a)","",
                            user.epsilon * user.secpera,&user.epsilon,&eps_set);CHKERRQ(ierr);
    if (eps_set) {  user.epsilon *= 1.0 / user.secpera;  }
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* Create machinery for parallel grid management (DA), nonlinear solver (SNES), 
     and Vecs for fields (solution, RHS).  Note default Mx=20 is 
     number of grid points.  Also degrees of freedom = 2 (thickness and velocity
     at each point) and stencil radius = ghost width = 1.                                    */
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,-20,2,1,PETSC_NULL,&user.da);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(user.da,0.0,user.xc,
                                 PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = DASetFieldName(user.da,0,"ice thickness (m)"); CHKERRQ(ierr);
  ierr = DASetFieldName(user.da,1,"ice velocity (m a-1)"); CHKERRQ(ierr);
  ierr = DAGetInfo(user.da,PETSC_IGNORE,&user.Mx,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,
                   PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
  ierr = DAGetCorners(user.da,&user.xs,PETSC_NULL,PETSC_NULL,&user.xm,PETSC_NULL,PETSC_NULL);
                   CHKERRQ(ierr);

  /* another DA for scalar parameters */
  ierr = DACreate1d(PETSC_COMM_WORLD,DA_NONPERIODIC,-20,1,1,PETSC_NULL,&user.scalarda);CHKERRQ(ierr);
  ierr = DASetUniformCoordinates(user.scalarda,0.0,user.xc,
                                 PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);

  /* Extract/allocate global vectors from DAs and duplicate for remaining same types */
  ierr = DACreateGlobalVector(user.da,&Hu);CHKERRQ(ierr);
  ierr = VecDuplicate(Hu,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(Hu,&user.Huexact);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.scalarda,&user.M);CHKERRQ(ierr);
  ierr = VecDuplicate(user.M,&user.B_stag);CHKERRQ(ierr);
  ierr = VecDuplicate(user.M,&user.beta);CHKERRQ(ierr);

  ierr = DASetLocalFunction(user.da,(DALocalFunction1)BodFunctionLocal);CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,r,SNESDAFormFunction,&user);CHKERRQ(ierr);

  /* setting up a matrix is only actually needed for -snes_fd case */
  ierr = DAGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,SNESDAComputeJacobian,PETSC_NULL);CHKERRQ(ierr);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  /* the the Bodvardsson (1955) exact solution allows setting M(x), B(x), beta(x), T(xc) */
  ierr = FillDistributedParams(&user);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"  using exact solution as initial guess\n");
             CHKERRQ(ierr);
  /* the exact thickness and exact ice velocity (user.uHexact) are known from Bodvardsson (1955) */
  ierr = FillExactSoln(&user); CHKERRQ(ierr);
  ierr = ScaleVelocity(&user,user.Huexact); CHKERRQ(ierr);
  ierr = VecCopy(user.Huexact,Hu); CHKERRQ(ierr);
  
  /************ SOLVE NONLINEAR SYSTEM  ************/
  /* recall that RHS  r  is used internally by KSP, and is set by the SNES */
  ierr = SNESSolve(snes,PETSC_NULL,Hu);CHKERRQ(ierr);

  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  %s Number of Newton iterations = %D\n",
           SNESConvergedReasons[reason],its);CHKERRQ(ierr);

  if (dump) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined result Hu\n");CHKERRQ(ierr);
    ierr = VecView(Hu,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing combined exact result Huexact\n");CHKERRQ(ierr);
    ierr = VecView(user.Huexact,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  viewing final combined residual at Hu\n");CHKERRQ(ierr);
    ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  }

  /* evaluate error relative to exact solution */
  ierr = MaxErrorNorms(&user,Hu,user.Huexact,&errHinf,&erruinf);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "  max norm errors:\n"
           "     |H - Hexact|_inf = %.4e m\n"
           "     |u - uexact|_inf = %.4e m a-1\n",
           errHinf,erruinf);CHKERRQ(ierr);

  ierr = VecDestroy(Hu);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = VecDestroy(user.Huexact);CHKERRQ(ierr);
  ierr = VecDestroy(user.M);CHKERRQ(ierr);
  ierr = VecDestroy(user.B_stag);CHKERRQ(ierr);
  ierr = VecDestroy(user.beta);CHKERRQ(ierr);

  ierr = MatDestroy(J); CHKERRQ(ierr);

  ierr = SNESDestroy(snes);CHKERRQ(ierr);

  ierr = DADestroy(user.da);CHKERRQ(ierr);
  ierr = DADestroy(user.scalarda);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

