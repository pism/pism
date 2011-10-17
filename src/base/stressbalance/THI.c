// Copyright (C) 2010-2011 Jed Brown and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

//! \file THI.c  Implementation of abstract data type (C language) for Jed's Blatter solver ("Toy Hydrostatic Ice").

#include "THItools.h"

static PetscCookie THI_COOKIE;

typedef enum {THIASSEMBLY_TRIDIAGONAL,THIASSEMBLY_FULL} THIAssemblyMode;

typedef struct _p_THI   *THI;

struct _p_THI {
  PETSCHEADER(int);
  void (*initialize)(THI,PetscReal x,PetscReal y,PrmNode *p);
  PetscInt  nlevels;
  PetscInt  zlevels;
  PetscReal Lx,Ly,Lz;           /* Model domain */
  PetscReal alpha;              /* Bed angle */
  Units     units;
  PetscReal dirichlet_scale;
  PetscReal ssa_friction_scale;
  PRange    eta;
  PRange    beta2;
  struct {
    PetscReal Bd2,eps,exponent;
  } viscosity;
  struct {
    PetscReal irefgam,eps2,exponent;
  } friction;
  PetscReal rhog;
  PetscTruth no_slip;
  PetscTruth tridiagonal;
  PetscTruth coarse2d;
  PetscTruth verbose;
  MatType mattype;
};

/* Tests A and C are from the ISMIP-HOM paper (Pattyn et al. 2008) */
static void THIInitialize_HOM_A(THI thi,PetscReal x,PetscReal y,PrmNode *p)
{
  Units units = thi->units;
  PetscReal s = -x*tan(thi->alpha);
  p->b = s - 1000*units->meter + 500*units->meter * sin(x*2*PETSC_PI/thi->Lx) * sin(y*2*PETSC_PI/thi->Ly);
  p->h = s - p->b;
  p->beta2 = 1e30;
}

static void THIInitialize_HOM_C(THI thi,PetscReal x,PetscReal y,PrmNode *p)
{
  Units units = thi->units;
  PetscReal s = -x*tan(thi->alpha);
  p->b = s - 1000*units->meter;
  p->h = s - p->b;
  /* tau_b = beta2 v   is a stress (Pa) */
  p->beta2 = 1000 * (1 + sin(x*2*PETSC_PI/thi->Lx)*sin(y*2*PETSC_PI/thi->Ly)) * units->Pascal * units->year / units->meter;
}

/* Same bed as test A, free slip everywhere except for a discontinuous jump to a circular sticky region in the middle. */
static void THIInitialize_HOM_X(THI thi,PetscReal xx,PetscReal yy,PrmNode *p)
{
  Units units = thi->units;
  PetscReal x = xx*2*PETSC_PI/thi->Lx - PETSC_PI,y = yy*2*PETSC_PI/thi->Ly - PETSC_PI; /* [-pi,pi] */
  PetscReal r = sqrt(x*x + y*y),s = -x*tan(thi->alpha);
  p->b = s - 1000*units->meter + 500*units->meter * sin(x + PETSC_PI) * sin(y + PETSC_PI);
  p->h = s - p->b;
  p->beta2 = 1000 * (r < 1 ? 2 : 0) * units->Pascal * units->year / units->meter;
}

/* Same bed as A, smoothly varying slipperiness, similar to Matlab's "sombrero" (uncorrelated with bathymetry) */
static void THIInitialize_HOM_Z(THI thi,PetscReal xx,PetscReal yy,PrmNode *p)
{
  Units units = thi->units;
  PetscReal x = xx*2*PETSC_PI/thi->Lx - PETSC_PI,y = yy*2*PETSC_PI/thi->Ly - PETSC_PI; /* [-pi,pi] */
  PetscReal r = sqrt(x*x + y*y),s = -x*tan(thi->alpha);
  p->b = s - 1000*units->meter + 500*units->meter * sin(x + PETSC_PI) * sin(y + PETSC_PI);
  p->h = s - p->b;
  p->beta2 = 1000 * (1. + sin(sqrt(16*r))/sqrt(1e-2 + 16*r)*cos(x*3/2)*cos(y*3/2)) * units->Pascal * units->year / units->meter;
}

static void THIFriction(THI thi,PetscReal rbeta2,PetscReal gam,PetscReal *beta2,PetscReal *dbeta2)
{
  if (thi->friction.irefgam == 0) {
    Units units = thi->units;
    thi->friction.irefgam = 1./(0.5*PetscSqr(100 * units->meter / units->year));
    thi->friction.eps2 = 0.5*PetscSqr(1.e-4 / thi->friction.irefgam);
  }
  if (thi->friction.exponent == 0) {
    *beta2 = rbeta2;
    *dbeta2 = 0;
  } else {
    *beta2 = rbeta2 * pow(thi->friction.eps2 + gam*thi->friction.irefgam,thi->friction.exponent);
    *dbeta2 = thi->friction.exponent * *beta2 / (thi->friction.eps2 + gam*thi->friction.irefgam) * thi->friction.irefgam;
  }
}

static void THIViscosity(THI thi,PetscReal gam,PetscReal *eta,PetscReal *deta)
{
  PetscReal Bd2,eps,exponent;
  if (thi->viscosity.Bd2 == 0) {
    Units units = thi->units;
    const PetscReal
      n = 3.,                                           /* Glen exponent */
      p = 1. + 1./n,                                    /* for Stokes */
      A = 1.e-16 * pow(units->Pascal,-n) / units->year, /* softness parameter (Pa^{-n}/s) */
      B = pow(A,-1./n);                                 /* hardness parameter */
    thi->viscosity.Bd2      = B/2;
    thi->viscosity.exponent = (p-2)/2;
    thi->viscosity.eps      = 0.5*PetscSqr(1e-5 / units->year);
  }
  Bd2      = thi->viscosity.Bd2;
  exponent = thi->viscosity.exponent;
  eps      = thi->viscosity.eps;
  *eta = Bd2 * pow(eps + gam,exponent);
  *deta = exponent * (*eta) / (eps + gam);
}


#undef __FUNCT__  
#define __FUNCT__ "THIDestroy"
static PetscErrorCode THIDestroy(THI thi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (--((PetscObject)thi)->refct > 0) PetscFunctionReturn(0);
  ierr = PetscFree(thi->units);CHKERRQ(ierr);
  ierr = PetscFree(thi->mattype);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(thi);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THICreate"
static PetscErrorCode THICreate(MPI_Comm comm,THI *inthi)
{
  static PetscTruth registered = PETSC_FALSE;
  THI thi;
  Units units;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *inthi = 0;
  if (!registered) {
    ierr = PetscCookieRegister("Toy Hydrostatic Ice",&THI_COOKIE);CHKERRQ(ierr);
    registered = PETSC_TRUE;
  }
  ierr = PetscHeaderCreate(thi,_p_THI,0,THI_COOKIE,-1,"THI",comm,THIDestroy,0);CHKERRQ(ierr);

  ierr = PetscNew(struct _n_Units,&thi->units);CHKERRQ(ierr);
  units = thi->units;
  units->meter  = 1e-2;
  units->second = 1e-7;
  units->kilogram = 1e-12;
  ierr = PetscOptionsBegin(comm,NULL,"Scaled units options","");CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-units_meter","1 meter in scaled length units","",units->meter,&units->meter,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-units_second","1 second in scaled time units","",units->second,&units->second,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-units_kilogram","1 kilogram in scaled mass units","",units->kilogram,&units->kilogram,NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  units->Pascal = units->kilogram / (units->meter * PetscSqr(units->second));
  units->year = 31556926. * units->second, /* seconds per year */

  thi->Lx              = 10.e3;
  thi->Ly              = 10.e3;
  thi->Lz              = 1000;
  thi->nlevels         = 1;
  thi->dirichlet_scale = 1;
  thi->verbose         = PETSC_FALSE;

  ierr = PetscOptionsBegin(comm,NULL,"Toy Hydrostatic Ice options","");CHKERRQ(ierr);
  {
    QuadratureType quad = QUAD_GAUSS;
    char homexp[] = "A";
    char mtype[256] = MATSBAIJ;
    PetscReal L,m = 1.0;
    PetscTruth flg;
    L = thi->Lx;
    ierr = PetscOptionsReal("-thi_L","Domain size (m)","",L,&L,&flg);CHKERRQ(ierr);
    if (flg) thi->Lx = thi->Ly = L;
    ierr = PetscOptionsReal("-thi_Lx","X Domain size (m)","",thi->Lx,&thi->Lx,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thi_Ly","Y Domain size (m)","",thi->Ly,&thi->Ly,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thi_Lz","Z Domain size (m)","",thi->Lz,&thi->Lz,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-thi_hom","ISMIP-HOM experiment (A or C)","",homexp,homexp,sizeof(homexp),NULL);CHKERRQ(ierr);
    switch (homexp[0] = (char)toupper(homexp[0])) {
      case 'A':
        thi->initialize = THIInitialize_HOM_A;
        thi->no_slip = PETSC_TRUE;
        thi->alpha = 0.5;
        break;
      case 'C':
        thi->initialize = THIInitialize_HOM_C;
        thi->no_slip = PETSC_FALSE;
        thi->alpha = 0.1;
        break;
      case 'X':
        thi->initialize = THIInitialize_HOM_X;
        thi->no_slip = PETSC_FALSE;
        thi->alpha = 0.3;
        break;
      case 'Z':
        thi->initialize = THIInitialize_HOM_Z;
        thi->no_slip = PETSC_FALSE;
        thi->alpha = 0.5;
        break;
      default:
        SETERRQ1(PETSC_ERR_SUP,"HOM experiment '%c' not implemented",homexp[0]);
    }
    ierr = PetscOptionsEnum("-thi_quadrature","Quadrature to use for 3D elements","",QuadratureTypes,(PetscEnum)quad,(PetscEnum*)&quad,NULL);CHKERRQ(ierr);
    switch (quad) {
      case QUAD_GAUSS:
        HexQInterp = HexQInterp_Gauss;
        HexQDeriv  = HexQDeriv_Gauss;
        break;
      case QUAD_LOBATTO:
        HexQInterp = HexQInterp_Lobatto;
        HexQDeriv  = HexQDeriv_Lobatto;
        break;
    }
    ierr = PetscOptionsReal("-thi_alpha","Bed angle (degrees)","",thi->alpha,&thi->alpha,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thi_friction_m","Friction exponent, 0=Coulomb, 1=Navier","",m,&m,NULL);CHKERRQ(ierr);
    thi->friction.exponent = (m-1)/2;
    ierr = PetscOptionsReal("-thi_dirichlet_scale","Scale Dirichlet boundary conditions by this factor","",thi->dirichlet_scale,&thi->dirichlet_scale,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thi_ssa_friction_scale","Scale slip boundary conditions by this factor in SSA (2D) assembly","",thi->ssa_friction_scale,&thi->ssa_friction_scale,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-thi_nlevels","Number of levels of refinement","",thi->nlevels,&thi->nlevels,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-thi_coarse2d","Use a 2D coarse space corresponding to SSA","",thi->coarse2d,&thi->coarse2d,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-thi_tridiagonal","Assemble a tridiagonal system (column coupling only) on the finest level","",thi->tridiagonal,&thi->tridiagonal,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsList("-thi_mat_type","Matrix type","MatSetType",MatList,mtype,(char*)mtype,sizeof(mtype),NULL);CHKERRQ(ierr);
    ierr = PetscStrallocpy(mtype,&thi->mattype);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-thi_verbose","Enable verbose output (like matrix sizes and statistics)","",thi->verbose,&thi->verbose,NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* dimensionalize */
  thi->Lx     *= units->meter;
  thi->Ly     *= units->meter;
  thi->Lz     *= units->meter;
  thi->alpha  *= PETSC_PI / 180;

  PRangeClear(&thi->eta);
  PRangeClear(&thi->beta2);

  {
    PetscReal u = 1000*units->meter/(3e7*units->second),
      gradu = u / (100*units->meter),eta,deta,
      rho = 910 * units->kilogram/pow(units->meter,3),
      grav = 9.81 * units->meter/PetscSqr(units->second),
      driving = rho * grav * tan(thi->alpha) * 1000*units->meter;
    THIViscosity(thi,0.5*gradu*gradu,&eta,&deta);
    thi->rhog = rho * grav;
    if (thi->verbose) {
      ierr = PetscPrintf(((PetscObject)thi)->comm,"Units: meter %8.2g  second %8.2g  kg %8.2g  Pa %8.2g\n",units->meter,units->second,units->kilogram,units->Pascal);CHKERRQ(ierr);
      ierr = PetscPrintf(((PetscObject)thi)->comm,"Domain (%6.2g,%6.2g,%6.2g), pressure %8.2g, driving stress %8.2g\n",thi->Lx,thi->Ly,thi->Lz,rho*grav*1e3*units->meter,driving);CHKERRQ(ierr);
      ierr = PetscPrintf(((PetscObject)thi)->comm,"Large velocity 1km/a %8.2g, velocity gradient %8.2g, eta %8.2g, stress %8.2g, ratio %8.2g\n",u,gradu,eta,2*eta*gradu,2*eta*gradu/driving);CHKERRQ(ierr);
      THIViscosity(thi,0.5*PetscSqr(1e-3*gradu),&eta,&deta);
      ierr = PetscPrintf(((PetscObject)thi)->comm,"Small velocity 1m/a  %8.2g, velocity gradient %8.2g, eta %8.2g, stress %8.2g, ratio %8.2g\n",1e-3*u,1e-3*gradu,eta,2*eta*1e-3*gradu,2*eta*1e-3*gradu/driving);CHKERRQ(ierr);
    }
  }

  *inthi = thi;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIInitializePrm"
static PetscErrorCode THIInitializePrm(THI thi,DA da2prm,Vec prm)
{
  PrmNode **p;
  PetscInt i,j,xs,xm,ys,ym,mx,my;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DAGetGhostCorners(da2prm,&ys,&xs,0,&ym,&xm,0);CHKERRQ(ierr);
  ierr = DAGetInfo(da2prm,0, &my,&mx,0, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
  ierr = DAVecGetArray(da2prm,prm,&p);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PetscReal xx = thi->Lx*i/mx,yy = thi->Ly*j/my;
      thi->initialize(thi,xx,yy,&p[i][j]);
    }
  }
  ierr = DAVecRestoreArray(da2prm,prm,&p);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THISetDMMG"
static PetscErrorCode THISetDMMG(THI thi,DMMG *dmmg)
{
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;
  if (DMMGGetLevels(dmmg) != thi->nlevels) SETERRQ(PETSC_ERR_ARG_CORRUPT,"DMMG nlevels does not agree with THI");
  for (i=0; i<thi->nlevels; i++) {
    PetscInt Mx,My,Mz,mx,my,s,dim;
    DAStencilType  st;
    DA da = (DA)dmmg[i]->dm,da2prm;
    Vec X;
    ierr = DAGetInfo(da,&dim, &Mz,&My,&Mx, 0,&my,&mx, 0,&s,0,&st);CHKERRQ(ierr);
    if (dim == 2) {
      ierr = DAGetInfo(da,&dim, &My,&Mx,0, &my,&mx,0, 0,&s,0,&st);CHKERRQ(ierr);
    }
    ierr = DACreate2d(((PetscObject)thi)->comm,DA_XYPERIODIC,st,My,Mx,my,mx,sizeof(PrmNode)/sizeof(PetscScalar),s,0,0,&da2prm);CHKERRQ(ierr);
    ierr = DACreateLocalVector(da2prm,&X);CHKERRQ(ierr);
    {
      PetscReal Lx = thi->Lx / thi->units->meter,Ly = thi->Ly / thi->units->meter,Lz = thi->Lz / thi->units->meter;
      if (dim == 2) {
        ierr = PetscPrintf(((PetscObject)thi)->comm,"Level %d domain size (m) %8.2g x %8.2g, num elements %3d x %3d (%8d), size (m) %g x %g\n",i,Lx,Ly,Mx,My,Mx*My,Lx/Mx,Ly/My);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(((PetscObject)thi)->comm,"Level %d domain size (m) %8.2g x %8.2g x %8.2g, num elements %3d x %3d x %3d (%8d), size (m) %g x %g x %g\n",i,Lx,Ly,Lz,Mx,My,Mz,Mx*My*Mz,Lx/Mx,Ly/My,1000./(Mz-1));CHKERRQ(ierr);
      }
    }
    ierr = THIInitializePrm(thi,da2prm,X);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)da,"DA2Prm",(PetscObject)da2prm);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)da,"DA2Prm_Vec",(PetscObject)X);CHKERRQ(ierr);
    ierr = DADestroy(da2prm);CHKERRQ(ierr);
    ierr = VecDestroy(X);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}


#undef __FUNCT__  
#define __FUNCT__ "THIInitial"
static PetscErrorCode THIInitial(DMMG dmmg,Vec X)
{
  THI         thi   = (THI)dmmg->user;
  DA          da    = (DA)dmmg->dm;
  PetscInt    i,j,k,xs,xm,ys,ym,zs,zm,mx,my;
  PetscReal   hx;
  PrmNode     **prm;
  Node        ***x;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DAGetInfo(da,0, 0,&my,&mx, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&zs,&ys,&xs,&zm,&ym,&xm);CHKERRQ(ierr);
  ierr = DAVecGetArray(da,X,&x);CHKERRQ(ierr);
  ierr = DAGetPrmNodeArray(da,&prm);CHKERRQ(ierr);
  hx = thi->Lx / mx;
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      for (k=zs; k<zs+zm; k++) {
        const PetscScalar zm1 = zm-1,
          drivingx = thi->rhog * (prm[i+1][j].b+prm[i+1][j].h - prm[i-1][j].b-prm[i-1][j].h) / (2*hx),
          drivingy = thi->rhog * (prm[i][j+1].b+prm[i][j+1].h - prm[i][j-1].b-prm[i][j-1].h) / (2*hx);
        x[i][j][k].u = 0. * drivingx * prm[i][j].h*(PetscScalar)k/zm1;
        x[i][j][k].v = 0. * drivingy * prm[i][j].h*(PetscScalar)k/zm1;
      }
    }
  }
  ierr = DAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  ierr = DARestorePrmNodeArray(da,&prm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static void PointwiseNonlinearity(THI thi,const Node n[restrict 8],const PetscReal phi[restrict 3],PetscReal dphi[restrict 8][3],PetscScalar *restrict u,PetscScalar *restrict v,PetscScalar du[restrict 3],PetscScalar dv[restrict 3],PetscReal *eta,PetscReal *deta)
{
  PetscInt l,ll;
  PetscScalar gam;

  du[0] = du[1] = du[2] = 0;
  dv[0] = dv[1] = dv[2] = 0;
  *u = 0;
  *v = 0;
  for (l=0; l<8; l++) {
    *u += phi[l] * n[l].u;
    *v += phi[l] * n[l].v;
    for (ll=0; ll<3; ll++) {
      du[ll] += dphi[l][ll] * n[l].u;
      dv[ll] += dphi[l][ll] * n[l].v;
    }
  }
  gam = Sqr(du[0]) + Sqr(dv[1]) + du[0]*dv[1] + 0.25*Sqr(du[1]+dv[0]) + 0.25*Sqr(du[2]) + 0.25*Sqr(dv[2]);
  THIViscosity(thi,PetscRealPart(gam),eta,deta);
}

static void PointwiseNonlinearity2D(THI thi,Node n[],PetscReal phi[],PetscReal dphi[4][2],PetscScalar *u,PetscScalar *v,PetscScalar du[],PetscScalar dv[],PetscReal *eta,PetscReal *deta)
{
  PetscInt l,ll;
  PetscScalar gam;

  du[0] = du[1] = 0;
  dv[0] = dv[1] = 0;
  *u = 0;
  *v = 0;
  for (l=0; l<4; l++) {
    *u += phi[l] * n[l].u;
    *v += phi[l] * n[l].v;
    for (ll=0; ll<2; ll++) {
      du[ll] += dphi[l][ll] * n[l].u;
      dv[ll] += dphi[l][ll] * n[l].v;
    }
  }
  gam = Sqr(du[0]) + Sqr(dv[1]) + du[0]*dv[1] + 0.25*Sqr(du[1]+dv[0]);
  THIViscosity(thi,PetscRealPart(gam),eta,deta);
}

#undef __FUNCT__  
#define __FUNCT__ "THIFunctionLocal"
static PetscErrorCode THIFunctionLocal(DALocalInfo *info,Node ***x,Node ***f,THI thi)
{
  PetscInt       xs,ys,xm,ym,zm,i,j,k,q,l;
  PetscReal      hx,hy,etamin,etamax,beta2min,beta2max;
  PrmNode        **prm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  xs = info->zs;
  ys = info->ys;
  xm = info->zm;
  ym = info->ym;
  zm = info->xm;
  hx = thi->Lx / info->mz;
  hy = thi->Ly / info->my;

  etamin   = 1e100;
  etamax   = 0;
  beta2min = 1e100;
  beta2max = 0;

  ierr = DAGetPrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PrmNode pn[4];
      QuadExtract(prm,i,j,pn);
      for (k=0; k<zm-1; k++) {
        PetscInt ls = 0;
        Node n[8],*fn[8];
        PetscReal zn[8],etabase = 0;
        PrmNodeHexGetZ(pn,k,zm,zn);
        HexExtract(x,i,j,k,n);
        HexExtractRef(f,i,j,k,fn);
        if (thi->no_slip && k == 0) {
          for (l=0; l<4; l++) n[l].u = n[l].v = 0;
          /* The first 4 basis functions lie on the bottom layer, so their contribution is exactly 0, hence we can skip them */
          ls = 4;
        }
        for (q=0; q<8; q++) {
          PetscReal dz[3],phi[8],dphi[8][3],jw,eta,deta;
          PetscScalar du[3],dv[3],u,v;
          HexGrad(HexQDeriv[q],zn,dz);
          HexComputeGeometry(q,hx,hy,dz,phi,dphi,&jw);
          PointwiseNonlinearity(thi,n,phi,dphi,&u,&v,du,dv,&eta,&deta);
          jw /= thi->rhog;      /* scales residuals to be O(1) */
          if (q == 0) etabase = eta;
          RangeUpdate(&etamin,&etamax,eta);
          for (l=ls; l<8; l++) { /* test functions */
            const PetscReal ds[2] = {-tan(thi->alpha),0};
            const PetscReal pp=phi[l],*dp = dphi[l];
            fn[l]->u += dp[0]*jw*eta*(4.*du[0]+2.*dv[1]) + dp[1]*jw*eta*(du[1]+dv[0]) + dp[2]*jw*eta*du[2] + pp*jw*thi->rhog*ds[0];
            fn[l]->v += dp[1]*jw*eta*(2.*du[0]+4.*dv[1]) + dp[0]*jw*eta*(du[1]+dv[0]) + dp[2]*jw*eta*dv[2] + pp*jw*thi->rhog*ds[1];
          }
        }
        if (k == 0) { /* we are on a bottom face */
          if (thi->no_slip) {
            /* Note: Non-Galerkin coarse grid operators are very sensitive to the scaling of Dirichlet boundary
            * conditions.  After shenanigans above, etabase contains the effective viscosity at the closest quadrature
            * point to the bed.  We want the diagonal entry in the Dirichlet condition to have similar magnitude to the
            * diagonal entry corresponding to the adjacent node.  The fundamental scaling of the viscous part is in
            * diagu, diagv below.  This scaling is easy to recognize by considering the finite difference operator after
            * scaling by element size.  The no-slip Dirichlet condition is scaled by this factor, and also in the
            * assembled matrix (see the similar block in THIJacobianLocal).
            */
            const PetscReal hz = PetscRealPart(pn[0].h)/(zm-1.);
            const PetscScalar diagu = 2*etabase/thi->rhog*(hx*hy/hz + hx*hz/hy + 4*hy*hz/hx),diagv = 2*etabase/thi->rhog*(hx*hy/hz + 4*hx*hz/hy + hy*hz/hx);
            fn[0]->u = thi->dirichlet_scale*diagu*n[0].u;
            fn[0]->v = thi->dirichlet_scale*diagv*n[0].v;
          } else {              /* Integrate over bottom face to apply boundary condition */
            for (q=0; q<4; q++) {
              const PetscReal jw = 0.25*hx*hy/thi->rhog,*phi = QuadQInterp[q];
              PetscScalar u=0,v=0,rbeta2=0;
              PetscReal beta2,dbeta2;
              for (l=0; l<4; l++) {
                u     += phi[l]*n[l].u;
                v     += phi[l]*n[l].v;
                rbeta2 += phi[l]*pn[l].beta2;
              }
              THIFriction(thi,PetscRealPart(rbeta2),PetscRealPart(u*u+v*v)/2,&beta2,&dbeta2);
              RangeUpdate(&beta2min,&beta2max,beta2);
              for (l=0; l<4; l++) {
                const PetscReal pp = phi[l];
                fn[ls+l]->u += pp*jw*beta2*u;
                fn[ls+l]->v += pp*jw*beta2*v;
              }
            }
          }
        }
      }
    }
  }

  ierr = DARestorePrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  ierr = PRangeMinMax(&thi->eta,etamin,etamax);CHKERRQ(ierr);
  ierr = PRangeMinMax(&thi->beta2,beta2min,beta2max);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIMatrixStatistics"
static PetscErrorCode THIMatrixStatistics(THI thi,Mat B,PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      nrm;
  PetscInt       m;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MatNorm(B,NORM_FROBENIUS,&nrm);CHKERRQ(ierr);
  ierr = MatGetSize(B,&m,0);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(((PetscObject)B)->comm,&rank);CHKERRQ(ierr);
  if (!rank) {
    PetscScalar val0,val2;
    ierr = MatGetValue(B,0,0,&val0);CHKERRQ(ierr);
    ierr = MatGetValue(B,2,2,&val2);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"Matrix dim %8d  norm %8.2e, (0,0) %8.2e  (2,2) %8.2e, eta [%8.2e,%8.2e] beta2 [%8.2e,%8.2e]\n",m,nrm,PetscRealPart(val0),PetscRealPart(val2),thi->eta.cmin,thi->eta.cmax,thi->beta2.cmin,thi->beta2.cmax);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THISurfaceStatistics"
static PetscErrorCode THISurfaceStatistics(DA da,Vec X,PetscReal *min,PetscReal *max,PetscReal *mean)
{
  PetscErrorCode ierr;
  Node           ***x;
  PetscInt       i,j,xs,ys,zs,xm,ym,zm,mx,my,mz;
  PetscReal      umin = 1e100,umax=-1e100;
  PetscScalar    usum=0.0,gusum;

  PetscFunctionBegin;
  *min = *max = *mean = 0;
  ierr = DAGetInfo(da,0, &mz,&my,&mx, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
  ierr = DAGetCorners(da,&zs,&ys,&xs,&zm,&ym,&xm);CHKERRQ(ierr);
  if (zs != 0 || zm != mz) SETERRQ(1,"Unexpected decomposition");
  ierr = DAVecGetArray(da,X,&x);CHKERRQ(ierr);
  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PetscReal u = PetscRealPart(x[i][j][zm-1].u);
      RangeUpdate(&umin,&umax,u);
      usum += u;
    }
  }
  ierr = DAVecRestoreArray(da,X,&x);CHKERRQ(ierr);
  ierr = PetscGlobalMin(&umin,min,((PetscObject)da)->comm);CHKERRQ(ierr);
  ierr = PetscGlobalMax(&umax,max,((PetscObject)da)->comm);CHKERRQ(ierr);
  ierr = PetscGlobalSum(&usum,&gusum,((PetscObject)da)->comm);CHKERRQ(ierr);
  *mean = PetscRealPart(gusum) / (mx*my);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THISolveStatistics"
static PetscErrorCode THISolveStatistics(THI thi,DMMG *dmmg,PetscInt coarsened,const char name[])
{
  MPI_Comm       comm    = ((PetscObject)thi)->comm;
  PetscInt       nlevels = DMMGGetLevels(dmmg),level = nlevels-1-coarsened;
  SNES           snes    = dmmg[level]->snes;
  Vec            X       = dmmg[level]->x;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscPrintf(comm,"Solution statistics after solve: %s\n",name);CHKERRQ(ierr);
  {
    PetscInt its,lits;
    SNESConvergedReason reason;
    ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);
    ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
    ierr = SNESGetLinearSolveIterations(snes,&lits);CHKERRQ(ierr);
    ierr = PetscPrintf(comm,"%s: Number of Newton iterations = %d, total linear iterations = %d\n",SNESConvergedReasons[reason],its,lits);CHKERRQ(ierr);
  }
  {
    PetscReal nrm2,min[3]={1e100,1e100,1e100},max[3]={-1e100,-1e100,-1e100};
    PetscInt i,j,m;
    PetscScalar *x;
    ierr = VecNorm(X,NORM_2,&nrm2);CHKERRQ(ierr);
    ierr = VecGetLocalSize(X,&m);CHKERRQ(ierr);
    ierr = VecGetArray(X,&x);CHKERRQ(ierr);
    for (i=0; i<m; i+=2) {
      PetscReal u = PetscRealPart(x[i]),v = PetscRealPart(x[i+1]),c = sqrt(u*u+v*v);
      min[0] = PetscMin(u,min[0]);
      min[1] = PetscMin(v,min[1]);
      min[2] = PetscMin(c,min[2]);
      max[0] = PetscMax(u,max[0]);
      max[1] = PetscMax(v,max[1]);
      max[2] = PetscMax(c,max[2]);
    }
    ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
    ierr = MPI_Allreduce(MPI_IN_PLACE,min,3,MPIU_REAL,MPI_MIN,((PetscObject)thi)->comm);CHKERRQ(ierr);
    ierr = MPI_Allreduce(MPI_IN_PLACE,max,3,MPIU_REAL,MPI_MAX,((PetscObject)thi)->comm);CHKERRQ(ierr);
    /* Dimensionalize to meters/year */
    nrm2 *= thi->units->year / thi->units->meter;
    for (j=0; j<3; j++) {
      min[j] *= thi->units->year / thi->units->meter;
      max[j] *= thi->units->year / thi->units->meter;
    }
    ierr = PetscPrintf(comm,"|X|_2 %g   u in [%g, %g]   v in [%g, %g]   c in [%g, %g] \n",nrm2,min[0],max[0],min[1],max[1],min[2],max[2]);CHKERRQ(ierr);
    {
      PetscReal umin,umax,umean;
      ierr = THISurfaceStatistics((DA)dmmg[level]->dm,X,&umin,&umax,&umean);CHKERRQ(ierr);
      umin  *= thi->units->year / thi->units->meter;
      umax  *= thi->units->year / thi->units->meter;
      umean *= thi->units->year / thi->units->meter;
      ierr = PetscPrintf(comm,"Surface statistics: u in [%12.6e, %12.6e] mean %12.6e\n",umin,umax,umean);CHKERRQ(ierr);
    }
    /* These values stay nondimensional */
    ierr = PetscPrintf(comm,"Global eta range   [%g, %g], converged range [%g, %g]\n",thi->eta.min,thi->eta.max,thi->eta.cmin,thi->eta.cmax);CHKERRQ(ierr);
    ierr = PetscPrintf(comm,"Global beta2 range [%g, %g], converged range [%g, %g]\n",thi->beta2.min,thi->beta2.max,thi->beta2.cmin,thi->beta2.cmax);CHKERRQ(ierr);
  }
  ierr = PetscPrintf(comm,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIJacobianLocal_2D"
static PetscErrorCode THIJacobianLocal_2D(DALocalInfo *info,Node **x,Mat B,THI thi)
{
  PetscInt       xs,ys,xm,ym,i,j,q,l,ll;
  PetscReal      hx,hy;
  PrmNode        **prm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  xs = info->ys;
  ys = info->xs;
  xm = info->ym;
  ym = info->xm;
  hx = thi->Lx / info->my;
  hy = thi->Ly / info->mx;

  ierr = MatZeroEntries(B);CHKERRQ(ierr);
  ierr = DAGetPrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      Node n[4];
      PrmNode pn[4];
      PetscScalar Ke[4*2][4*2];
      QuadExtract(prm,i,j,pn);
      QuadExtract(x,i,j,n);
      PetscMemzero(Ke,sizeof(Ke));
      for (q=0; q<4; q++) {
        PetscReal phi[4],dphi[4][2],jw,eta,deta,beta2,dbeta2;
        PetscScalar u,v,du[2],dv[2],h = 0,rbeta2 = 0;
        for (l=0; l<4; l++) {
          phi[l] = QuadQInterp[q][l];
          dphi[l][0] = QuadQDeriv[q][l][0]*2./hx;
          dphi[l][1] = QuadQDeriv[q][l][1]*2./hy;
          h += phi[l] * pn[l].h;
          rbeta2 += phi[l] * pn[l].beta2;
        }
        jw = 0.25*hx*hy / thi->rhog; /* rhog is only scaling */
        PointwiseNonlinearity2D(thi,n,phi,dphi,&u,&v,du,dv,&eta,&deta);
        THIFriction(thi,PetscRealPart(rbeta2),PetscRealPart(u*u+v*v)/2,&beta2,&dbeta2);
        for (l=0; l<4; l++) {
          const PetscReal pp = phi[l],*dp = dphi[l];
          for (ll=0; ll<4; ll++) {
            const PetscReal ppl = phi[ll],*dpl = dphi[ll];
            PetscScalar dgdu,dgdv;
            dgdu = 2.*du[0]*dpl[0] + dv[1]*dpl[0] + 0.5*(du[1]+dv[0])*dpl[1];
            dgdv = 2.*dv[1]*dpl[1] + du[0]*dpl[1] + 0.5*(du[1]+dv[0])*dpl[0];
            /* Picard part */
            Ke[l*2+0][ll*2+0] += dp[0]*jw*eta*4.*dpl[0] + dp[1]*jw*eta*dpl[1] + pp*jw*(beta2/h)*ppl*thi->ssa_friction_scale;
            Ke[l*2+0][ll*2+1] += dp[0]*jw*eta*2.*dpl[1] + dp[1]*jw*eta*dpl[0];
            Ke[l*2+1][ll*2+0] += dp[1]*jw*eta*2.*dpl[0] + dp[0]*jw*eta*dpl[1];
            Ke[l*2+1][ll*2+1] += dp[1]*jw*eta*4.*dpl[1] + dp[0]*jw*eta*dpl[0] + pp*jw*(beta2/h)*ppl*thi->ssa_friction_scale;
            /* extra Newton terms */
            Ke[l*2+0][ll*2+0] += dp[0]*jw*deta*dgdu*(4.*du[0]+2.*dv[1]) + dp[1]*jw*deta*dgdu*(du[1]+dv[0]) + pp*jw*(dbeta2/h)*u*u*ppl*thi->ssa_friction_scale;
            Ke[l*2+0][ll*2+1] += dp[0]*jw*deta*dgdv*(4.*du[0]+2.*dv[1]) + dp[1]*jw*deta*dgdv*(du[1]+dv[0]) + pp*jw*(dbeta2/h)*u*v*ppl*thi->ssa_friction_scale;
            Ke[l*2+1][ll*2+0] += dp[1]*jw*deta*dgdu*(4.*dv[1]+2.*du[0]) + dp[0]*jw*deta*dgdu*(du[1]+dv[0]) + pp*jw*(dbeta2/h)*v*u*ppl*thi->ssa_friction_scale;
            Ke[l*2+1][ll*2+1] += dp[1]*jw*deta*dgdv*(4.*dv[1]+2.*du[0]) + dp[0]*jw*deta*dgdv*(du[1]+dv[0]) + pp*jw*(dbeta2/h)*v*v*ppl*thi->ssa_friction_scale;
          }
        }
      }
      {
        const MatStencil rc[4] = {{0,i,j,0},{0,i+1,j,0},{0,i+1,j+1,0},{0,i,j+1,0}};
        ierr = MatSetValuesBlockedStencil(B,4,rc,4,rc,&Ke[0][0],ADD_VALUES);CHKERRQ(ierr);
      }
    }
  }
  ierr = DARestorePrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(B,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  if (thi->verbose) {ierr = THIMatrixStatistics(thi,B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIJacobianLocal_3D"
static PetscErrorCode THIJacobianLocal_3D(DALocalInfo *info,Node ***x,Mat B,THI thi,THIAssemblyMode amode)
{
  PetscInt       xs,ys,xm,ym,zm,i,j,k,q,l,ll;
  PetscReal      hx,hy;
  PrmNode        **prm;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  xs = info->zs;
  ys = info->ys;
  xm = info->zm;
  ym = info->ym;
  zm = info->xm;
  hx = thi->Lx / info->mz;
  hy = thi->Ly / info->my;

  ierr = MatZeroEntries(B);CHKERRQ(ierr);
  ierr = DAGetPrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PrmNode pn[4];
      QuadExtract(prm,i,j,pn);
      for (k=0; k<zm-1; k++) {
        Node n[8];
        PetscReal zn[8],etabase = 0;
        PetscScalar Ke[8*2][8*2];
        PetscInt ls = 0;

        PrmNodeHexGetZ(pn,k,zm,zn);
        HexExtract(x,i,j,k,n);
        PetscMemzero(Ke,sizeof(Ke));
        if (thi->no_slip && k == 0) {
          for (l=0; l<4; l++) n[l].u = n[l].v = 0;
          ls = 4;
        }
        for (q=0; q<8; q++) {
          PetscReal dz[3],phi[8],dphi[8][3],jw,eta,deta;
          PetscScalar du[3],dv[3],u,v;
          HexGrad(HexQDeriv[q],zn,dz);
          HexComputeGeometry(q,hx,hy,dz,phi,dphi,&jw);
          PointwiseNonlinearity(thi,n,phi,dphi,&u,&v,du,dv,&eta,&deta);
          jw /= thi->rhog;      /* residuals are scaled by this factor */
          if (q == 0) etabase = eta;
          for (l=ls; l<8; l++) { /* test functions */
            const PetscReal *restrict dp = dphi[l];
#if USE_SSE2_KERNELS
            /* gcc (up to my 4.5 snapshot) is really bad at hoisting intrinsics so we do it manually */
            __m128d
              p4 = _mm_set1_pd(4),p2 = _mm_set1_pd(2),p05 = _mm_set1_pd(0.5),
              p42 = _mm_setr_pd(4,2),p24 = _mm_shuffle_pd(p42,p42,_MM_SHUFFLE2(0,1)),
              du0 = _mm_set1_pd(du[0]),du1 = _mm_set1_pd(du[1]),du2 = _mm_set1_pd(du[2]),
              dv0 = _mm_set1_pd(dv[0]),dv1 = _mm_set1_pd(dv[1]),dv2 = _mm_set1_pd(dv[2]),
              jweta = _mm_set1_pd(jw*eta),jwdeta = _mm_set1_pd(jw*deta),
              dp0 = _mm_set1_pd(dp[0]),dp1 = _mm_set1_pd(dp[1]),dp2 = _mm_set1_pd(dp[2]),
              dp0jweta = _mm_mul_pd(dp0,jweta),dp1jweta = _mm_mul_pd(dp1,jweta),dp2jweta = _mm_mul_pd(dp2,jweta),
              p4du0p2dv1 = _mm_add_pd(_mm_mul_pd(p4,du0),_mm_mul_pd(p2,dv1)), /* 4 du0 + 2 dv1 */
              p4dv1p2du0 = _mm_add_pd(_mm_mul_pd(p4,dv1),_mm_mul_pd(p2,du0)), /* 4 dv1 + 2 du0 */
              pdu2dv2 = _mm_unpacklo_pd(du2,dv2),                             /* [du2, dv2] */
              du1pdv0 = _mm_add_pd(du1,dv0),                                  /* du1 + dv0 */
              t1 = _mm_mul_pd(dp0,p4du0p2dv1),                                /* dp0 (4 du0 + 2 dv1) */
              t2 = _mm_mul_pd(dp1,p4dv1p2du0);                                /* dp1 (4 dv1 + 2 du0) */

#endif
#if defined COMPUTE_LOWER_TRIANGULAR  /* The element matrices are always symmetric so computing the lower-triangular part is not necessary */
            for (ll=ls; ll<8; ll++) { /* trial functions */
#else
            for (ll=l; ll<8; ll++) {
#endif
              const PetscReal *restrict dpl = dphi[ll];
              if (amode == THIASSEMBLY_TRIDIAGONAL && (l-ll)%4) continue; /* these entries would not be inserted */
#if !USE_SSE2_KERNELS
              /* The analytic Jacobian in nice, easy-to-read form */
              {
                PetscScalar dgdu,dgdv;
                dgdu = 2.*du[0]*dpl[0] + dv[1]*dpl[0] + 0.5*(du[1]+dv[0])*dpl[1] + 0.5*du[2]*dpl[2];
                dgdv = 2.*dv[1]*dpl[1] + du[0]*dpl[1] + 0.5*(du[1]+dv[0])*dpl[0] + 0.5*dv[2]*dpl[2];
                /* Picard part */
                Ke[l*2+0][ll*2+0] += dp[0]*jw*eta*4.*dpl[0] + dp[1]*jw*eta*dpl[1] + dp[2]*jw*eta*dpl[2];
                Ke[l*2+0][ll*2+1] += dp[0]*jw*eta*2.*dpl[1] + dp[1]*jw*eta*dpl[0];
                Ke[l*2+1][ll*2+0] += dp[1]*jw*eta*2.*dpl[0] + dp[0]*jw*eta*dpl[1];
                Ke[l*2+1][ll*2+1] += dp[1]*jw*eta*4.*dpl[1] + dp[0]*jw*eta*dpl[0] + dp[2]*jw*eta*dpl[2];
                /* extra Newton terms */
                Ke[l*2+0][ll*2+0] += dp[0]*jw*deta*dgdu*(4.*du[0]+2.*dv[1]) + dp[1]*jw*deta*dgdu*(du[1]+dv[0]) + dp[2]*jw*deta*dgdu*du[2];
                Ke[l*2+0][ll*2+1] += dp[0]*jw*deta*dgdv*(4.*du[0]+2.*dv[1]) + dp[1]*jw*deta*dgdv*(du[1]+dv[0]) + dp[2]*jw*deta*dgdv*du[2];
                Ke[l*2+1][ll*2+0] += dp[1]*jw*deta*dgdu*(4.*dv[1]+2.*du[0]) + dp[0]*jw*deta*dgdu*(du[1]+dv[0]) + dp[2]*jw*deta*dgdu*dv[2];
                Ke[l*2+1][ll*2+1] += dp[1]*jw*deta*dgdv*(4.*dv[1]+2.*du[0]) + dp[0]*jw*deta*dgdv*(du[1]+dv[0]) + dp[2]*jw*deta*dgdv*dv[2];
              }
#else
              /* This SSE2 code is an exact replica of above, but uses explicit packed instructions for some speed
              * benefit.  On my hardware, these intrinsics are almost twice as fast as above, reducing total assembly cost
              * by 25 to 30 percent. */
              {
                __m128d
                  keu = _mm_loadu_pd(&Ke[l*2+0][ll*2+0]),
                  kev = _mm_loadu_pd(&Ke[l*2+1][ll*2+0]),
                  dpl01 = _mm_loadu_pd(&dpl[0]),dpl10 = _mm_shuffle_pd(dpl01,dpl01,_MM_SHUFFLE2(0,1)),dpl2 = _mm_set_sd(dpl[2]),
                  t0,t3,pdgduv;
                keu = _mm_add_pd(keu,_mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp0jweta,p42),dpl01),
                                                _mm_add_pd(_mm_mul_pd(dp1jweta,dpl10),
                                                           _mm_mul_pd(dp2jweta,dpl2))));
                kev = _mm_add_pd(kev,_mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp1jweta,p24),dpl01),
                                                _mm_add_pd(_mm_mul_pd(dp0jweta,dpl10),
                                                           _mm_mul_pd(dp2jweta,_mm_shuffle_pd(dpl2,dpl2,_MM_SHUFFLE2(0,1))))));
                pdgduv = _mm_mul_pd(p05,_mm_add_pd(_mm_add_pd(_mm_mul_pd(p42,_mm_mul_pd(du0,dpl01)),
                                                              _mm_mul_pd(p24,_mm_mul_pd(dv1,dpl01))),
                                                   _mm_add_pd(_mm_mul_pd(du1pdv0,dpl10),
                                                              _mm_mul_pd(pdu2dv2,_mm_set1_pd(dpl[2]))))); /* [dgdu, dgdv] */
                t0 = _mm_mul_pd(jwdeta,pdgduv);  /* jw deta [dgdu, dgdv] */
                t3 = _mm_mul_pd(t0,du1pdv0);     /* t0 (du1 + dv0) */
                _mm_storeu_pd(&Ke[l*2+0][ll*2+0],_mm_add_pd(keu,_mm_add_pd(_mm_mul_pd(t1,t0),
                                                                          _mm_add_pd(_mm_mul_pd(dp1,t3),
                                                                                     _mm_mul_pd(t0,_mm_mul_pd(dp2,du2))))));
                _mm_storeu_pd(&Ke[l*2+1][ll*2+0],_mm_add_pd(kev,_mm_add_pd(_mm_mul_pd(t2,t0),
                                                                          _mm_add_pd(_mm_mul_pd(dp0,t3),
                                                                                     _mm_mul_pd(t0,_mm_mul_pd(dp2,dv2))))));
              }
#endif
            }
          }
        }
        if (k == 0) { /* on a bottom face */
          if (thi->no_slip) {
            const PetscReal hz = PetscRealPart(pn[0].h)/(zm-1);
            const PetscScalar diagu = 2*etabase/thi->rhog*(hx*hy/hz + hx*hz/hy + 4*hy*hz/hx),diagv = 2*etabase/thi->rhog*(hx*hy/hz + 4*hx*hz/hy + hy*hz/hx);
            Ke[0][0] = thi->dirichlet_scale*diagu;
            Ke[1][1] = thi->dirichlet_scale*diagv;
          } else {
            for (q=0; q<4; q++) {
              const PetscReal jw = 0.25*hx*hy/thi->rhog,*phi = QuadQInterp[q];
              PetscScalar u=0,v=0,rbeta2=0;
              PetscReal beta2,dbeta2;
              for (l=0; l<4; l++) {
                u     += phi[l]*n[l].u;
                v     += phi[l]*n[l].v;
                rbeta2 += phi[l]*pn[l].beta2;
              }
              THIFriction(thi,PetscRealPart(rbeta2),PetscRealPart(u*u+v*v)/2,&beta2,&dbeta2);
              for (l=0; l<4; l++) {
                const PetscReal pp = phi[l];
                for (ll=0; ll<4; ll++) {
                  const PetscReal ppl = phi[ll];
                  Ke[l*2+0][ll*2+0] += pp*jw*beta2*ppl + pp*jw*dbeta2*u*u*ppl;
                  Ke[l*2+0][ll*2+1] +=                   pp*jw*dbeta2*u*v*ppl;
                  Ke[l*2+1][ll*2+0] +=                   pp*jw*dbeta2*v*u*ppl;
                  Ke[l*2+1][ll*2+1] += pp*jw*beta2*ppl + pp*jw*dbeta2*v*v*ppl;
                }
              }
            }
          }
        }
        {
          const MatStencil rc[8] = {{i,j,k,0},{i+1,j,k,0},{i+1,j+1,k,0},{i,j+1,k,0},{i,j,k+1,0},{i+1,j,k+1,0},{i+1,j+1,k+1,0},{i,j+1,k+1,0}};
          if (amode == THIASSEMBLY_TRIDIAGONAL) {
            for (l=0; l<4; l++) { /* Copy out each of the blocks, discarding horizontal coupling */
              const PetscInt l4 = l+4;
              const MatStencil rcl[2] = {{rc[l].k,rc[l].j,rc[l].i,0},{rc[l4].k,rc[l4].j,rc[l4].i,0}};
#if defined COMPUTE_LOWER_TRIANGULAR
              const PetscScalar Kel[4][4] = {{Ke[2*l+0][2*l+0] ,Ke[2*l+0][2*l+1] ,Ke[2*l+0][2*l4+0] ,Ke[2*l+0][2*l4+1]},
                                             {Ke[2*l+1][2*l+0] ,Ke[2*l+1][2*l+1] ,Ke[2*l+1][2*l4+0] ,Ke[2*l+1][2*l4+1]},
                                             {Ke[2*l4+0][2*l+0],Ke[2*l4+0][2*l+1],Ke[2*l4+0][2*l4+0],Ke[2*l4+0][2*l4+1]},
                                             {Ke[2*l4+1][2*l+0],Ke[2*l4+1][2*l+1],Ke[2*l4+1][2*l4+0],Ke[2*l4+1][2*l4+1]}};
#else
              /* Same as above except for the lower-left block */
              const PetscScalar Kel[4][4] = {{Ke[2*l+0][2*l+0] ,Ke[2*l+0][2*l+1] ,Ke[2*l+0][2*l4+0] ,Ke[2*l+0][2*l4+1]},
                                             {Ke[2*l+1][2*l+0] ,Ke[2*l+1][2*l+1] ,Ke[2*l+1][2*l4+0] ,Ke[2*l+1][2*l4+1]},
                                             {Ke[2*l+0][2*l4+0],Ke[2*l+1][2*l4+0],Ke[2*l4+0][2*l4+0],Ke[2*l4+0][2*l4+1]},
                                             {Ke[2*l+0][2*l4+1],Ke[2*l+1][2*l4+1],Ke[2*l4+1][2*l4+0],Ke[2*l4+1][2*l4+1]}};
#endif
              ierr = MatSetValuesBlockedStencil(B,2,rcl,2,rcl,&Kel[0][0],ADD_VALUES);CHKERRQ(ierr);
            }
          } else {
#if !defined COMPUTE_LOWER_TRIANGULAR /* fill in lower-triangular part, this is really cheap compared to computing the entries */
            for (l=0; l<8; l++) {
              for (ll=l+1; ll<8; ll++) {
                Ke[ll*2+0][l*2+0] = Ke[l*2+0][ll*2+0];
                Ke[ll*2+1][l*2+0] = Ke[l*2+0][ll*2+1];
                Ke[ll*2+0][l*2+1] = Ke[l*2+1][ll*2+0];
                Ke[ll*2+1][l*2+1] = Ke[l*2+1][ll*2+1];
              }
            }
#endif
            ierr = MatSetValuesBlockedStencil(B,8,rc,8,rc,&Ke[0][0],ADD_VALUES);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = DARestorePrmNodeArray(info->da,&prm);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(B,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  if (thi->verbose) {ierr = THIMatrixStatistics(thi,B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIJacobianLocal_3D_Full"
static PetscErrorCode THIJacobianLocal_3D_Full(DALocalInfo *info,Node ***x,Mat B,THI thi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = THIJacobianLocal_3D(info,x,B,thi,THIASSEMBLY_FULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIJacobianLocal_3D_Tridiagonal"
static PetscErrorCode THIJacobianLocal_3D_Tridiagonal(DALocalInfo *info,Node ***x,Mat B,THI thi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = THIJacobianLocal_3D(info,x,B,thi,THIASSEMBLY_TRIDIAGONAL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DARefineHierarchy_THI"
static PetscErrorCode DARefineHierarchy_THI(DA dac0,PetscInt nlevels,DA hierarchy[])
{
  PetscErrorCode ierr;
  THI thi;
  PetscInt dim,M,N,m,n,s,dof;
  DA dac,daf;
  DAStencilType  st;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)dac0,"THI",(PetscObject*)&thi);CHKERRQ(ierr);
  if (!thi) SETERRQ(PETSC_ERR_ARG_WRONG,"Cannot refine this DA, missing composed THI instance");
  if (nlevels > 1) {
    ierr = DARefineHierarchy(dac0,nlevels-1,hierarchy);CHKERRQ(ierr);
    dac = hierarchy[nlevels-2];
  } else {
    dac = dac0;
  }
  ierr = DAGetInfo(dac,&dim, &N,&M,0, &n,&m,0, &dof,&s,0,&st);CHKERRQ(ierr);
  if (dim != 2) SETERRQ(PETSC_ERR_ARG_WRONG,"This function can only refine 2D DAs");
  /* Creates a 3D DA with the same map-plane layout as the 2D one, with contiguous columns */
  ierr = DACreate3d(((PetscObject)dac)->comm,DA_YZPERIODIC,st,thi->zlevels,N,M,1,n,m,dof,s,NULL,NULL,NULL,&daf);CHKERRQ(ierr);
  daf->ops->getmatrix        = dac->ops->getmatrix;
  daf->ops->getinterpolation = dac->ops->getinterpolation;
  daf->ops->getcoloring      = dac->ops->getcoloring;
  daf->interptype            = dac->interptype;

  ierr = DASetFieldName(daf,0,"x-velocity");CHKERRQ(ierr);
  ierr = DASetFieldName(daf,1,"y-velocity");CHKERRQ(ierr);
  hierarchy[nlevels-1] = daf;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAGetInterpolation_THI"
static PetscErrorCode DAGetInterpolation_THI(DA dac,DA daf,Mat *A,Vec *scale)
{
  PetscErrorCode ierr;
  PetscTruth flg,isda2;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dac,DM_COOKIE,1);
  PetscValidHeaderSpecific(daf,DM_COOKIE,2);
  PetscValidPointer(A,3);
  if (scale) {
    (void) ierr;   /* avoid compiler warning "empty body in an if-statement" */
    PetscValidPointer(scale,4);
  }
  ierr = PetscTypeCompare((PetscObject)dac,DA2D,&flg);
  if (!flg) SETERRQ(PETSC_ERR_ARG_WRONG,"Expected coarse DA to be 2D");
  ierr = PetscTypeCompare((PetscObject)daf,DA2D,&isda2);CHKERRQ(ierr);
  if (isda2) {
    /* We are in the 2D problem and use normal DA interpolation */
    ierr = DAGetInterpolation(dac,daf,A,scale);CHKERRQ(ierr);
  } else {
    PetscInt i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz,rstart,cstart;
    Mat B;

    ierr = DAGetInfo(daf,0, &mz,&my,&mx, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
    ierr = DAGetCorners(daf,&zs,&ys,&xs,&zm,&ym,&xm);CHKERRQ(ierr);
    if (zs != 0) SETERRQ(1,"unexpected");
    ierr = MatCreate(((PetscObject)daf)->comm,&B);CHKERRQ(ierr);
    ierr = MatSetSizes(B,xm*ym*zm,xm*ym,mx*my*mz,mx*my);CHKERRQ(ierr);
    
    ierr = MatSetType(B,MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(B,1,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(B,1,NULL,0,NULL);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(B,&rstart,NULL);CHKERRQ(ierr);
    ierr = MatGetOwnershipRangeColumn(B,&cstart,NULL);CHKERRQ(ierr);
    for (i=xs; i<xs+xm; i++) {
      for (j=ys; j<ys+ym; j++) {
        for (k=zs; k<zs+zm; k++) {
          PetscInt i2 = i*ym+j,i3 = i2*zm+k;
          PetscScalar val = ((k == 0 || k == mz-1) ? 0.5 : 1.) / (mz-1.); /* Integration using trapezoid rule */
          ierr = MatSetValue(B,cstart+i3,rstart+i2,val,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatCreateMAIJ(B,sizeof(Node)/sizeof(PetscScalar),A);CHKERRQ(ierr);
    ierr = MatDestroy(B);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "DAGetMatrix_THI_Tridiagonal"
static PetscErrorCode DAGetMatrix_THI_Tridiagonal(DA da,const MatType mtype,Mat *J)
{
  PetscErrorCode ierr;
  Mat A;
  PetscInt xm,ym,zm,dim,dof = 2,starts[3],dims[3];
  ISLocalToGlobalMapping ltog,ltogb;

  PetscFunctionBegin;
  ierr = DAGetInfo(da,&dim, 0,0,0, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
  if (dim != 3) SETERRQ(PETSC_ERR_ARG_WRONG,"Expected DA to be 3D");
  ierr = DAGetCorners(da,0,0,0,&zm,&ym,&xm);CHKERRQ(ierr);
  ierr = DAGetISLocalToGlobalMapping(da,&ltog);CHKERRQ(ierr);
  ierr = DAGetISLocalToGlobalMappingBlck(da,&ltogb);CHKERRQ(ierr);
  ierr = MatCreate(((PetscObject)da)->comm,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,dof*xm*ym*zm,dof*xm*ym*zm,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetType(A,mtype);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,6,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,6,PETSC_NULL,0,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSeqBAIJSetPreallocation(A,dof,3,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIBAIJSetPreallocation(A,dof,3,PETSC_NULL,0,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSeqSBAIJSetPreallocation(A,dof,2,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPISBAIJSetPreallocation(A,dof,2,PETSC_NULL,0,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSetBlockSize(A,dof);CHKERRQ(ierr);
  ierr = MatSetLocalToGlobalMapping(A,ltog);CHKERRQ(ierr);
  ierr = MatSetLocalToGlobalMappingBlock(A,ltogb);CHKERRQ(ierr);
  ierr = DAGetGhostCorners(da,&starts[0],&starts[1],&starts[2],&dims[0],&dims[1],&dims[2]);CHKERRQ(ierr);
  ierr = MatSetStencil(A,dim,dims,starts,dof);CHKERRQ(ierr);
  *J = A;
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "THIDAVecView_VTK_XML"
static PetscErrorCode THIDAVecView_VTK_XML(THI thi,DA da,Vec X,const char filename[])
{
  const PetscInt dof   = 2;
  Units          units = thi->units;
  MPI_Comm       comm;
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscMPIInt    rank,size,tag,nn,nmax;
  PetscInt       mx,my,mz,r,range[6];
  PetscScalar    *x;

  PetscFunctionBegin;
  comm = ((PetscObject)thi)->comm;
  ierr = DAGetInfo(da,0, &mz,&my,&mx, 0,0,0, 0,0,0,0);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(comm,filename,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",0,mz-1,0,my-1,0,mx-1);CHKERRQ(ierr);

  ierr = DAGetCorners(da,range,range+1,range+2,range+3,range+4,range+5);CHKERRQ(ierr);
  nn = PetscMPIIntCast(range[3]*range[4]*range[5]*dof);
  ierr = MPI_Reduce(&nn,&nmax,1,MPI_INT,MPI_MAX,0,comm);CHKERRQ(ierr);
  tag  = ((PetscObject) viewer)->tag;
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  if (!rank) {
    PetscScalar *array;
    ierr = PetscMalloc(nmax*sizeof(PetscScalar),&array);CHKERRQ(ierr);
    for (r=0; r<size; r++) {
      PetscInt i,j,k,xs,xm,ys,ym,zs,zm;
      PetscScalar *ptr;
      MPI_Status status;
      if (r) {
        ierr = MPI_Recv(range,6,MPIU_INT,r,tag,comm,MPI_STATUS_IGNORE);CHKERRQ(ierr);
      }
      zs = range[0];ys = range[1];xs = range[2];zm = range[3];ym = range[4];xm = range[5];
      if (xm*ym*zm*dof > nmax) SETERRQ(1,"should not happen");
      if (r) {
        ierr = MPI_Recv(array,nmax,MPIU_SCALAR,r,tag,comm,&status);CHKERRQ(ierr);
        ierr = MPI_Get_count(&status,MPIU_SCALAR,&nn);CHKERRQ(ierr);
        if (nn != xm*ym*zm*dof) SETERRQ(1,"should not happen");
        ptr = array;
      } else ptr = x;
      ierr = PetscViewerASCIIPrintf(viewer,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",zs,zs+zm-1,ys,ys+ym-1,xs,xs+xm-1);CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(viewer,"      <Points>\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");CHKERRQ(ierr);
      for (i=xs; i<xs+xm; i++) {
        for (j=ys; j<ys+ym; j++) {
          for (k=zs; k<zs+zm; k++) {
            PrmNode p;
            PetscReal xx = thi->Lx*i/mx,yy = thi->Ly*j/my,zz;
            thi->initialize(thi,xx,yy,&p);
            zz = PetscRealPart(p.b) + PetscRealPart(p.h)*k/(mz-1);
            ierr = PetscViewerASCIIPrintf(viewer,"%f %f %f\n",xx,yy,zz);CHKERRQ(ierr);
          }
        }
      }
      ierr = PetscViewerASCIIPrintf(viewer,"        </DataArray>\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"      </Points>\n");CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(viewer,"      <PointData>\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");CHKERRQ(ierr);
      for (i=0; i<nn; i+=dof) {
        ierr = PetscViewerASCIIPrintf(viewer,"%f %f %f\n",PetscRealPart(ptr[i])*units->year/units->meter,PetscRealPart(ptr[i+1])*units->year/units->meter,0.0);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"        </DataArray>\n");CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(viewer,"        <DataArray type=\"Int32\" Name=\"rank\" NumberOfComponents=\"1\" format=\"ascii\">\n");CHKERRQ(ierr);
      for (i=0; i<nn; i+=dof) {
        ierr = PetscViewerASCIIPrintf(viewer,"%d\n",r);CHKERRQ(ierr);
      }
      ierr = PetscViewerASCIIPrintf(viewer,"        </DataArray>\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"      </PointData>\n");CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(viewer,"    </Piece>\n");CHKERRQ(ierr);
    }
    ierr = PetscFree(array);CHKERRQ(ierr);
  } else {
    ierr = MPI_Send(range,6,MPIU_INT,0,tag,comm);CHKERRQ(ierr);
    ierr = MPI_Send(x,nn,MPIU_SCALAR,0,tag,comm);CHKERRQ(ierr);
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"  </StructuredGrid>\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"</VTKFile>\n");CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

