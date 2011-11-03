/* Copyright (C) 2010-2011 Jed Brown, Ed Bueler and Constantine Khroulev */

/* This file is part of PISM. */

/* PISM is free software; you can redistribute it and/or modify it under the */
/* terms of the GNU General Public License as published by the Free Software */
/* Foundation; either version 2 of the License, or (at your option) any later */
/* version. */

/* PISM is distributed in the hope that it will be useful, but WITHOUT ANY */
/* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS */
/* FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more */
/* details. */

/* You should have received a copy of the GNU General Public License */
/* along with PISM; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */

/* ! \file THI.c  Implementation of abstract data type (C language) for Jed's Blatter solver ("Toy Hydrostatic Ice"). */

#include "THI.hh"
#include "THItools.hh"
static PetscCookie THI_COOKIE;

struct _p_THI {
  PETSCHEADER(int);
  PetscInt  nlevels;            // number of multigrid levels
  PetscInt  zlevels;            // number of vertical levels
  PetscReal Lx, Ly;             /* Model domain */
  PetscReal alpha;              /* Bed angle */
  Units     units;
  PetscReal dirichlet_scale;
  PRange    eta;
  PRange    beta2;
  struct {
    PetscReal Bd2, eps, exponent;
  } viscosity;
  struct {
    PetscReal irefgam, eps2, exponent;
  } friction;
  PetscReal rhog;
  PetscTruth no_slip;
  PetscTruth tridiagonal;
  PetscTruth verbose;
  MatType mattype;
};

//! \brief Implements equation for the friction parameter beta^2.
/*!
 * \f[ \beta^2(\gamma_b) = \beta^2_0 ( \epsilon_b^2 + \gamma_b )^{\frac{m-1}{2}} \f]
 *
 * beta2(gamma_b) = beta2_0 (eps_b^2 + gamma_b)^((1-m)/2)
 *
 * All the inputs (rbeta2, gam) are evaluated at quadrature points.
 * gam (gamma_b) is 0.5 * (u^2 + v^2).
 * rbeta2 is the spatially-variable friction parameter corresponding to beta2_0.
 *
 * Note that in PISM (eq. 18 in [\ref BBssasliding], regulatized):
 *
 * \f[ \tau_{b,i} = -\tau_c \frac{v_i}{|\mathbf{v}|} = -\tau_c v_i (\epsilon + |\mathbf{v}|^2)^(-1/2). \f]
 *
 * Compare to \f$\beta^2\f$ in [\ref BrownSmithAhmadia] with \f$m=0\f$:
 *
 * \f[\beta^2(\gamma_b) = \beta_0 \left[\frac12 (\epsilon_b^2 + |\mathbf{v}|^2)\right]^{-1/2}  =
 * \beta_0\sqrt{2}(\epsilon_b^2 + |\mathbf{v}|^2)^{-1/2}. \f]
 */
static void THIFriction(THI thi, PetscReal rbeta2, PetscReal gam, PetscReal *beta2, PetscReal *dbeta2)
{
  // Compute THI friction parameters during the first call:
  if (thi->friction.irefgam == 0) {
    Units units = thi->units;
    thi->friction.irefgam = 1./(0.5*PetscSqr(100 * units->meter / units->year));
    thi->friction.eps2 = 0.5*PetscSqr(1.e-4 / thi->friction.irefgam);
  }
  // Use cached irefgam and eps2:
  if (thi->friction.exponent == 0) {
    *beta2 = rbeta2;
    *dbeta2 = 0;
  } else {
    *beta2 = rbeta2 * pow(thi->friction.eps2 + gam*thi->friction.irefgam, thi->friction.exponent);
    *dbeta2 = thi->friction.exponent * *beta2 / (thi->friction.eps2 + gam*thi->friction.irefgam) * thi->friction.irefgam;
  }
}

//! \brief Implements equation (2.2) in [\ref BrownSmithAhmadia].
/*!
 * \f[ \eta(\gamma) = \frac{B}{2}\left(\epsilon^2/2 + \gamma\right)^{\frac{1-n}{2n}} \f]
 *
 * eta(gamma) = B/2 * (eps^2/s + gamma)^((1-n)/(2n))
 *
 * FIXME: We need to allow spatially-variable ice hardness B.
 */
static void THIViscosity(THI thi, PetscReal gam, PetscReal *eta, PetscReal *deta)
{
  PetscReal Bd2, eps, exponent;
  // Compute THI viscosity parameters during the first call:
  if (thi->viscosity.Bd2 == 0) {
    Units units = thi->units;
    const PetscReal
      n = 3., /* Glen exponent */
      p = 1. + 1./n, /* for Stokes */
      A = 1.e-16 * pow(units->Pascal, -n) / units->year, /* softness parameter (Pa^{-n}/s) */
      B = pow(A, -1./n);                                 /* hardness parameter */
    thi->viscosity.Bd2      = B/2;
    thi->viscosity.exponent = (p-2)/2;
    thi->viscosity.eps      = 0.5*PetscSqr(1e-5 / units->year);
  }

  // Use cached Bd2, exponent, eps:
  Bd2      = thi->viscosity.Bd2;
  exponent = thi->viscosity.exponent;
  eps      = thi->viscosity.eps;
  *eta = Bd2 * pow(eps + gam, exponent);
  *deta = exponent * (*eta) / (eps + gam);
}


#undef __FUNCT__
#define __FUNCT__ "THIDestroy"
//! \brief De-allocate the THI object.
PetscErrorCode THIDestroy(THI thi)
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
//! \brief Allocate the THI object.
PetscErrorCode THICreate(MPI_Comm comm, THI *inthi)
{
  static PetscTruth registered = PETSC_FALSE;
  THI thi;
  Units units;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *inthi = 0;
  if (!registered) {
    ierr = PetscCookieRegister("Toy Hydrostatic Ice", &THI_COOKIE);CHKERRQ(ierr);
    registered = PETSC_TRUE;
  }
  ierr = PetscHeaderCreate(thi, _p_THI, 0, THI_COOKIE, -1, "THI", comm, THIDestroy, 0);CHKERRQ(ierr);

  ierr = PetscNew(struct _n_Units, &thi->units);CHKERRQ(ierr);
  units = thi->units;
  units->meter  = 1e-2;
  units->second = 1e-7;
  units->kilogram = 1e-12;
  ierr = PetscOptionsBegin(comm, NULL, "Scaled units options", "");CHKERRQ(ierr);
  {
    ierr = PetscOptionsReal("-units_meter", "1 meter in scaled length units", "", units->meter, &units->meter, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-units_second", "1 second in scaled time units", "", units->second, &units->second, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-units_kilogram", "1 kilogram in scaled mass units", "", units->kilogram, &units->kilogram, NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  units->Pascal = units->kilogram / (units->meter * PetscSqr(units->second));
  units->year = 31556926. * units->second, /* seconds per year */

  thi->nlevels         = 1;
  thi->dirichlet_scale = 1;
  thi->verbose         = PETSC_FALSE;

  thi->alpha = 0;
  thi->no_slip = PETSC_TRUE;

  ierr = PetscOptionsBegin(comm, NULL, "Toy Hydrostatic Ice options", "");CHKERRQ(ierr);
  {
    QuadratureType quad = QUAD_GAUSS;
    char mtype[256] = MATSBAIJ;
    PetscReal m = 1.0;

    ierr = PetscOptionsEnum("-thi_quadrature", "Quadrature to use for 3D elements", "", QuadratureTypes, (PetscEnum)quad, (PetscEnum*)&quad, NULL);CHKERRQ(ierr);
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
    ierr = PetscOptionsReal("-thi_alpha", "Bed angle (degrees)", "", thi->alpha, &thi->alpha, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thi_friction_m", "Friction exponent, 0=Coulomb, 1=Navier", "", m, &m, NULL);CHKERRQ(ierr);
    thi->friction.exponent = (m-1)/2;
    ierr = PetscOptionsReal("-thi_dirichlet_scale", "Scale Dirichlet boundary conditions by this factor", "", thi->dirichlet_scale, &thi->dirichlet_scale, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-thi_nlevels", "Number of levels of refinement", "", thi->nlevels, &thi->nlevels, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-thi_tridiagonal", "Assemble a tridiagonal system (column coupling only) on the finest level", "", thi->tridiagonal, &thi->tridiagonal, NULL);CHKERRQ(ierr);
    ierr = PetscOptionsList("-thi_mat_type", "Matrix type", "MatSetType", MatList, mtype, (char*)mtype, sizeof(mtype), NULL);CHKERRQ(ierr);
    ierr = PetscStrallocpy(mtype, &thi->mattype);CHKERRQ(ierr);
    ierr = PetscOptionsTruth("-thi_verbose", "Enable verbose output (like matrix sizes and statistics)", "", thi->verbose, &thi->verbose, NULL);CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  /* dimensionalize */
  thi->alpha  *= PETSC_PI / 180;

  PRangeClear(&thi->eta);
  PRangeClear(&thi->beta2);

  {
    PetscReal u = 1000*units->meter/(3e7*units->second),
      gradu = u / (100*units->meter), eta, deta,
      rho = 910 * units->kilogram/pow(units->meter, 3),
      grav = 9.81 * units->meter/PetscSqr(units->second),
      driving = rho * grav * tan(thi->alpha) * 1000*units->meter;
    THIViscosity(thi, 0.5*gradu*gradu, &eta, &deta);
    thi->rhog = rho * grav;
    if (thi->verbose) {
      ierr = PetscPrintf(((PetscObject)thi)->comm, "Units: meter %8.2g  second %8.2g  kg %8.2g  Pa %8.2g\n", units->meter, units->second, units->kilogram, units->Pascal);CHKERRQ(ierr);
      ierr = PetscPrintf(((PetscObject)thi)->comm, "Large velocity 1km/a %8.2g, velocity gradient %8.2g, eta %8.2g, stress %8.2g, ratio %8.2g\n", u, gradu, eta, 2*eta*gradu, 2*eta*gradu/driving);CHKERRQ(ierr);
      THIViscosity(thi, 0.5*PetscSqr(1e-3*gradu), &eta, &deta);
      ierr = PetscPrintf(((PetscObject)thi)->comm, "Small velocity 1m/a  %8.2g, velocity gradient %8.2g, eta %8.2g, stress %8.2g, ratio %8.2g\n", 1e-3*u, 1e-3*gradu, eta, 2*eta*1e-3*gradu, 2*eta*1e-3*gradu/driving);CHKERRQ(ierr);
    }
  }

  *inthi = thi;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "THISetDMMG"
PetscErrorCode THISetDMMG(THI thi, DMMG *dmmg)
{
  PetscErrorCode ierr;
  PetscInt i;

  PetscFunctionBegin;
  if (DMMGGetLevels(dmmg) != thi->nlevels) SETERRQ(PETSC_ERR_ARG_CORRUPT, "DMMG nlevels does not agree with THI");
  for (i=0; i<thi->nlevels; i++) {
    PetscInt Mx, My, Mz, mx, my, s, dim;
    DAStencilType  st;
    DA da = (DA)dmmg[i]->dm, da2prm;
    Vec X;
    ierr = DAGetInfo(da, &dim, &Mz, &My, &Mx, 0, &my, &mx, 0, &s, 0, &st);CHKERRQ(ierr);
    ierr = DACreate2d(((PetscObject)thi)->comm,
                      DA_XYPERIODIC,
                      st, // stencil type
                      My, Mx, // number of grid points in each direction
                      my, mx, // number of processors in each direction
                      sizeof(PrmNode)/sizeof(PetscScalar), // number of degrees of freedom
                      s, // stencil width
                      PETSC_NULL, PETSC_NULL, // pointers to lists defining partitions in each direction
                      &da2prm);CHKERRQ(ierr);
    ierr = DACreateLocalVector(da2prm, &X);CHKERRQ(ierr);
    {
      PetscReal Lx = thi->Lx / thi->units->meter, Ly = thi->Ly / thi->units->meter;
      ierr = PetscPrintf(((PetscObject)thi)->comm,
                         "Level %d domain size (m) %8.2g x %8.2g,"
                         " num elements %3d x %3d x %3d (%8d), size (m) %g x %g\n",
                         i, Lx, Ly, Mx, My, Mz, Mx*My*Mz, Lx/Mx, Ly/My);CHKERRQ(ierr);
    }
    ierr = PetscObjectCompose((PetscObject)da, "DA2Prm", (PetscObject)da2prm);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)da, "DA2Prm_Vec", (PetscObject)X);CHKERRQ(ierr);
    ierr = DADestroy(da2prm);CHKERRQ(ierr);
    ierr = VecDestroy(X);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static void PointwiseNonlinearity(THI thi, const Node n[restrict 8], const PetscReal phi[restrict 3], PetscReal dphi[restrict 8][3], PetscScalar *restrict u, PetscScalar *restrict v, PetscScalar du[restrict 3], PetscScalar dv[restrict 3], PetscReal *eta, PetscReal *deta)
{
  PetscInt l, ll;
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
  THIViscosity(thi, PetscRealPart(gam), eta, deta);
}

#undef __FUNCT__
#define __FUNCT__ "THIFunctionLocal"
PetscErrorCode THIFunctionLocal(DALocalInfo *info, Node ***x, Node ***f, THI thi)
{
  PetscInt       xs, ys, xm, ym, zm, i, j, k, q, l;
  PetscReal      hx, hy, etamin, etamax, beta2min, beta2max;
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

  ierr = DAGetPrmNodeArray(info->da, &prm);CHKERRQ(ierr);

  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PrmNode pn[4];
      QuadExtract(prm, i, j, pn);
      for (k=0; k<zm-1; k++) {
        PetscInt ls = 0;
        Node n[8], *fn[8];
        PetscReal zn[8], etabase = 0;
        PrmNodeHexGetZ(pn, k, zm, zn);
        HexExtract(x, i, j, k, n);
        HexExtractRef(f, i, j, k, fn);
        if (thi->no_slip && k == 0) {
          for (l=0; l<4; l++) n[l].u = n[l].v = 0;
          /* The first 4 basis functions lie on the bottom layer, so their contribution is exactly 0, hence we can skip them */
          ls = 4;
        }
        for (q=0; q<8; q++) {
          PetscReal dz[3], phi[8], dphi[8][3], jw, eta, deta;
          PetscScalar du[3], dv[3], u, v;
          HexGrad(HexQDeriv[q], zn, dz);
          HexComputeGeometry(q, hx, hy, dz, phi, dphi, &jw);
          PointwiseNonlinearity(thi, n, phi, dphi, &u, &v, du, dv, &eta, &deta);
          jw /= thi->rhog;      /* scales residuals to be O(1) */
          if (q == 0) etabase = eta;
          RangeUpdate(&etamin, &etamax, eta);
          for (l=ls; l<8; l++) { /* test functions */
            const PetscReal ds[2] = {-tan(thi->alpha), 0};
            const PetscReal pp=phi[l], *dp = dphi[l];
            fn[l]->u += dp[0]*jw*eta*(4.*du[0]+2.*dv[1]) + dp[1]*jw*eta*(du[1]+dv[0]) + dp[2]*jw*eta*du[2] + pp*jw*thi->rhog*ds[0];
            fn[l]->v += dp[1]*jw*eta*(2.*du[0]+4.*dv[1]) + dp[0]*jw*eta*(du[1]+dv[0]) + dp[2]*jw*eta*dv[2] + pp*jw*thi->rhog*ds[1];
          }
        }
        if (k == 0) { /* we are on a bottom face */
          if (thi->no_slip) {
            /* Note: Non-Galerkin coarse grid operators are very sensitive to
            * the scaling of Dirichlet boundary conditions. After shenanigans
            * above, etabase contains the effective viscosity at the closest
            * quadrature point to the bed. We want the diagonal entry in the
            * Dirichlet condition to have similar magnitude to the diagonal
            * entry corresponding to the adjacent node. The fundamental scaling
            * of the viscous part is in diagu, diagv below. This scaling is
            * easy to recognize by considering the finite difference operator
            * after scaling by element size. The no-slip Dirichlet condition is
            * scaled by this factor, and also in the assembled matrix (see the
            * similar block in THIJacobianLocal).
            */
            const PetscReal hz = PetscRealPart(pn[0].h)/(zm-1.);
            const PetscScalar diagu = 2*etabase/thi->rhog*(hx*hy/hz + hx*hz/hy + 4*hy*hz/hx),
              diagv = 2*etabase/thi->rhog*(hx*hy/hz + 4*hx*hz/hy + hy*hz/hx);
            fn[0]->u = thi->dirichlet_scale*diagu*n[0].u;
            fn[0]->v = thi->dirichlet_scale*diagv*n[0].v;
          } else {              /* Integrate over bottom face to apply boundary condition */
            for (q=0; q<4; q++) {
              const PetscReal jw = 0.25*hx*hy/thi->rhog, *phi = QuadQInterp[q];
              PetscScalar u=0, v=0, rbeta2=0;
              PetscReal beta2, dbeta2;
              for (l=0; l<4; l++) {
                u     += phi[l]*n[l].u;
                v     += phi[l]*n[l].v;
                rbeta2 += phi[l]*pn[l].beta2;
              }
              THIFriction(thi, PetscRealPart(rbeta2), PetscRealPart(u*u+v*v)/2, &beta2, &dbeta2);
              RangeUpdate(&beta2min, &beta2max, beta2);
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

  ierr = DARestorePrmNodeArray(info->da, &prm);CHKERRQ(ierr);

  ierr = PRangeMinMax(&thi->eta, etamin, etamax);CHKERRQ(ierr);
  ierr = PRangeMinMax(&thi->beta2, beta2min, beta2max);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "THIMatrixStatistics"
static PetscErrorCode THIMatrixStatistics(THI thi, Mat B, PetscViewer viewer)
{
  PetscErrorCode ierr;
  PetscReal      nrm;
  PetscInt       m;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MatNorm(B, NORM_FROBENIUS, &nrm);CHKERRQ(ierr);
  ierr = MatGetSize(B, &m, 0);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(((PetscObject)B)->comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    PetscScalar val0, val2;
    ierr = MatGetValue(B, 0, 0, &val0);CHKERRQ(ierr);
    ierr = MatGetValue(B, 2, 2, &val2);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, "Matrix dim %8d  norm %8.2e, (0, 0) %8.2e  (2, 2) %8.2e, eta [%8.2e, %8.2e] beta2 [%8.2e, %8.2e]\n", m, nrm, PetscRealPart(val0), PetscRealPart(val2), thi->eta.cmin, thi->eta.cmax, thi->beta2.cmin, thi->beta2.cmax);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "THIJacobianLocal_3D"
PetscErrorCode THIJacobianLocal_3D(DALocalInfo *info, Node ***x, Mat B, THI thi, THIAssemblyMode amode)
{
  PetscInt       xs, ys, xm, ym, zm, i, j, k, q, l, ll;
  PetscReal      hx, hy;
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
  ierr = DAGetPrmNodeArray(info->da, &prm);CHKERRQ(ierr);

  for (i=xs; i<xs+xm; i++) {
    for (j=ys; j<ys+ym; j++) {
      PrmNode pn[4];
      QuadExtract(prm, i, j, pn);
      for (k=0; k<zm-1; k++) {
        Node n[8];
        PetscReal zn[8], etabase = 0;
        PetscScalar Ke[8*2][8*2];
        PetscInt ls = 0;

        PrmNodeHexGetZ(pn, k, zm, zn);
        HexExtract(x, i, j, k, n);
        PetscMemzero(Ke, sizeof(Ke));
        if (thi->no_slip && k == 0) {
          for (l=0; l<4; l++) n[l].u = n[l].v = 0;
          ls = 4;
        }
        for (q=0; q<8; q++) {
          PetscReal dz[3], phi[8], dphi[8][3], jw, eta, deta;
          PetscScalar du[3], dv[3], u, v;
          HexGrad(HexQDeriv[q], zn, dz);
          HexComputeGeometry(q, hx, hy, dz, phi, dphi, &jw);
          PointwiseNonlinearity(thi, n, phi, dphi, &u, &v, du, dv, &eta, &deta);
          jw /= thi->rhog;      /* residuals are scaled by this factor */
          if (q == 0) etabase = eta;
          for (l=ls; l<8; l++) { /* test functions */
            const PetscReal *restrict dp = dphi[l];
#if USE_SSE2_KERNELS
            /* gcc (up to my 4.5 snapshot) is really bad at hoisting intrinsics so we do it manually */
            __m128d
              p4 = _mm_set1_pd(4), p2 = _mm_set1_pd(2), p05 = _mm_set1_pd(0.5),
              p42 = _mm_setr_pd(4, 2), p24 = _mm_shuffle_pd(p42, p42, _MM_SHUFFLE2(0, 1)),
              du0 = _mm_set1_pd(du[0]), du1 = _mm_set1_pd(du[1]), du2 = _mm_set1_pd(du[2]),
              dv0 = _mm_set1_pd(dv[0]), dv1 = _mm_set1_pd(dv[1]), dv2 = _mm_set1_pd(dv[2]),
              jweta = _mm_set1_pd(jw*eta), jwdeta = _mm_set1_pd(jw*deta),
              dp0 = _mm_set1_pd(dp[0]), dp1 = _mm_set1_pd(dp[1]), dp2 = _mm_set1_pd(dp[2]),
              dp0jweta = _mm_mul_pd(dp0, jweta), dp1jweta = _mm_mul_pd(dp1, jweta), dp2jweta = _mm_mul_pd(dp2, jweta),
              p4du0p2dv1 = _mm_add_pd(_mm_mul_pd(p4, du0), _mm_mul_pd(p2, dv1)), /* 4 du0 + 2 dv1 */
              p4dv1p2du0 = _mm_add_pd(_mm_mul_pd(p4, dv1), _mm_mul_pd(p2, du0)), /* 4 dv1 + 2 du0 */
              pdu2dv2 = _mm_unpacklo_pd(du2, dv2), /* [du2, dv2] */
              du1pdv0 = _mm_add_pd(du1, dv0), /* du1 + dv0 */
              t1 = _mm_mul_pd(dp0, p4du0p2dv1), /* dp0 (4 du0 + 2 dv1) */
              t2 = _mm_mul_pd(dp1, p4dv1p2du0);                                /* dp1 (4 dv1 + 2 du0) */

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
                PetscScalar dgdu, dgdv;
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
                  dpl01 = _mm_loadu_pd(&dpl[0]), dpl10 = _mm_shuffle_pd(dpl01, dpl01, _MM_SHUFFLE2(0, 1)), dpl2 = _mm_set_sd(dpl[2]),
                  t0, t3, pdgduv;
                keu = _mm_add_pd(keu, _mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp0jweta, p42), dpl01),
                                                _mm_add_pd(_mm_mul_pd(dp1jweta, dpl10),
                                                           _mm_mul_pd(dp2jweta, dpl2))));
                kev = _mm_add_pd(kev, _mm_add_pd(_mm_mul_pd(_mm_mul_pd(dp1jweta, p24), dpl01),
                                                _mm_add_pd(_mm_mul_pd(dp0jweta, dpl10),
                                                           _mm_mul_pd(dp2jweta, _mm_shuffle_pd(dpl2, dpl2, _MM_SHUFFLE2(0, 1))))));
                pdgduv = _mm_mul_pd(p05, _mm_add_pd(_mm_add_pd(_mm_mul_pd(p42, _mm_mul_pd(du0, dpl01)),
                                                              _mm_mul_pd(p24, _mm_mul_pd(dv1, dpl01))),
                                                   _mm_add_pd(_mm_mul_pd(du1pdv0, dpl10),
                                                              _mm_mul_pd(pdu2dv2, _mm_set1_pd(dpl[2]))))); /* [dgdu, dgdv] */
                t0 = _mm_mul_pd(jwdeta, pdgduv);  /* jw deta [dgdu, dgdv] */
                t3 = _mm_mul_pd(t0, du1pdv0);     /* t0 (du1 + dv0) */
                _mm_storeu_pd(&Ke[l*2+0][ll*2+0], _mm_add_pd(keu, _mm_add_pd(_mm_mul_pd(t1, t0),
                                                                          _mm_add_pd(_mm_mul_pd(dp1, t3),
                                                                                     _mm_mul_pd(t0, _mm_mul_pd(dp2, du2))))));
                _mm_storeu_pd(&Ke[l*2+1][ll*2+0], _mm_add_pd(kev, _mm_add_pd(_mm_mul_pd(t2, t0),
                                                                          _mm_add_pd(_mm_mul_pd(dp0, t3),
                                                                                     _mm_mul_pd(t0, _mm_mul_pd(dp2, dv2))))));
              }
#endif
            }
          }
        }
        if (k == 0) { /* on a bottom face */
          if (thi->no_slip) {
            const PetscReal hz = PetscRealPart(pn[0].h)/(zm-1);
            const PetscScalar diagu = 2*etabase/thi->rhog*(hx*hy/hz + hx*hz/hy + 4*hy*hz/hx), diagv = 2*etabase/thi->rhog*(hx*hy/hz + 4*hx*hz/hy + hy*hz/hx);
            Ke[0][0] = thi->dirichlet_scale*diagu;
            Ke[1][1] = thi->dirichlet_scale*diagv;
          } else {
            for (q=0; q<4; q++) {
              const PetscReal jw = 0.25*hx*hy/thi->rhog, *phi = QuadQInterp[q];
              PetscScalar u=0, v=0, rbeta2=0;
              PetscReal beta2, dbeta2;
              for (l=0; l<4; l++) {
                u     += phi[l]*n[l].u;
                v     += phi[l]*n[l].v;
                rbeta2 += phi[l]*pn[l].beta2;
              }
              THIFriction(thi, PetscRealPart(rbeta2), PetscRealPart(u*u+v*v)/2, &beta2, &dbeta2);
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
          const MatStencil rc[8] = {{i, j, k, 0}, {i+1, j, k, 0}, {i+1, j+1, k, 0}, {i, j+1, k, 0}, {i, j, k+1, 0}, {i+1, j, k+1, 0}, {i+1, j+1, k+1, 0}, {i, j+1, k+1, 0}};
          if (amode == THIASSEMBLY_TRIDIAGONAL) {
            for (l=0; l<4; l++) { /* Copy out each of the blocks, discarding horizontal coupling */
              const PetscInt l4 = l+4;
              const MatStencil rcl[2] = {{rc[l].k, rc[l].j, rc[l].i, 0}, {rc[l4].k, rc[l4].j, rc[l4].i, 0}};
#if defined COMPUTE_LOWER_TRIANGULAR
              const PetscScalar Kel[4][4] = {{Ke[2*l+0][2*l+0] , Ke[2*l+0][2*l+1] , Ke[2*l+0][2*l4+0] , Ke[2*l+0][2*l4+1]},
                                             {Ke[2*l+1][2*l+0] , Ke[2*l+1][2*l+1] , Ke[2*l+1][2*l4+0] , Ke[2*l+1][2*l4+1]},
                                             {Ke[2*l4+0][2*l+0], Ke[2*l4+0][2*l+1], Ke[2*l4+0][2*l4+0], Ke[2*l4+0][2*l4+1]},
                                             {Ke[2*l4+1][2*l+0], Ke[2*l4+1][2*l+1], Ke[2*l4+1][2*l4+0], Ke[2*l4+1][2*l4+1]}};
#else
              /* Same as above except for the lower-left block */
              const PetscScalar Kel[4][4] = {{Ke[2*l+0][2*l+0] , Ke[2*l+0][2*l+1] , Ke[2*l+0][2*l4+0] , Ke[2*l+0][2*l4+1]},
                                             {Ke[2*l+1][2*l+0] , Ke[2*l+1][2*l+1] , Ke[2*l+1][2*l4+0] , Ke[2*l+1][2*l4+1]},
                                             {Ke[2*l+0][2*l4+0], Ke[2*l+1][2*l4+0], Ke[2*l4+0][2*l4+0], Ke[2*l4+0][2*l4+1]},
                                             {Ke[2*l+0][2*l4+1], Ke[2*l+1][2*l4+1], Ke[2*l4+1][2*l4+0], Ke[2*l4+1][2*l4+1]}};
#endif
              ierr = MatSetValuesBlockedStencil(B, 2, rcl, 2, rcl, &Kel[0][0], ADD_VALUES);CHKERRQ(ierr);
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
            ierr = MatSetValuesBlockedStencil(B, 8, rc, 8, rc, &Ke[0][0], ADD_VALUES);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = DARestorePrmNodeArray(info->da, &prm);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetOption(B, MAT_SYMMETRIC, PETSC_TRUE);CHKERRQ(ierr);
  if (thi->verbose) {ierr = THIMatrixStatistics(thi, B, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "THIJacobianLocal_3D_Full"
PetscErrorCode THIJacobianLocal_3D_Full(DALocalInfo *info, Node ***x, Mat B, THI thi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = THIJacobianLocal_3D(info, x, B, thi, THIASSEMBLY_FULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "THIJacobianLocal_3D_Tridiagonal"
PetscErrorCode THIJacobianLocal_3D_Tridiagonal(DALocalInfo *info, Node ***x, Mat B, THI thi)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = THIJacobianLocal_3D(info, x, B, thi, THIASSEMBLY_TRIDIAGONAL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DARefineHierarchy_THI"
PetscErrorCode DARefineHierarchy_THI(DA dac0, PetscInt nlevels, DA hierarchy[])
{
  PetscErrorCode ierr;
  THI thi;
  PetscInt dim, M, N, m, n, s, dof;
  DA dac, daf;
  DAStencilType  st;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)dac0, "THI", (PetscObject*)&thi);CHKERRQ(ierr);
  if (!thi) SETERRQ(PETSC_ERR_ARG_WRONG, "Cannot refine this DA, missing composed THI instance");
  if (nlevels > 1) {
    ierr = DARefineHierarchy(dac0, nlevels-1, hierarchy);CHKERRQ(ierr);
    dac = hierarchy[nlevels-2];
  } else {
    dac = dac0;
  }
  ierr = DAGetInfo(dac, &dim, &N, &M, 0, &n, &m, 0, &dof, &s, 0, &st);CHKERRQ(ierr);
  if (dim != 2) SETERRQ(PETSC_ERR_ARG_WRONG, "This function can only refine 2D DAs");
  /* Creates a 3D DA with the same map-plane layout as the 2D one, with contiguous columns */
  ierr = DACreate3d(((PetscObject)dac)->comm, DA_YZPERIODIC, st, thi->zlevels, N, M, 1, n, m, dof, s, NULL, NULL, NULL, &daf);CHKERRQ(ierr);
  daf->ops->getmatrix        = dac->ops->getmatrix;
  daf->ops->getinterpolation = dac->ops->getinterpolation;
  daf->ops->getcoloring      = dac->ops->getcoloring;
  daf->interptype            = dac->interptype;

  ierr = DASetFieldName(daf, 0, "x-velocity");CHKERRQ(ierr);
  ierr = DASetFieldName(daf, 1, "y-velocity");CHKERRQ(ierr);
  hierarchy[nlevels-1] = daf;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAGetInterpolation_THI"
PetscErrorCode DAGetInterpolation_THI(DA dac, DA daf, Mat *A, Vec *scale)
{
  PetscErrorCode ierr;
  PetscTruth flg, isda2;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dac, DM_COOKIE, 1);
  PetscValidHeaderSpecific(daf, DM_COOKIE, 2);
  PetscValidPointer(A, 3);
  if (scale) {
    (void) ierr;   /* avoid compiler warning "empty body in an if-statement" */
    PetscValidPointer(scale, 4);
  }
  ierr = PetscTypeCompare((PetscObject)dac, DA2D, &flg);
  if (!flg) SETERRQ(PETSC_ERR_ARG_WRONG, "Expected coarse DA to be 2D");
  ierr = PetscTypeCompare((PetscObject)daf, DA2D, &isda2);CHKERRQ(ierr);
  if (isda2) {
    /* We are in the 2D problem and use normal DA interpolation */
    ierr = DAGetInterpolation(dac, daf, A, scale);CHKERRQ(ierr);
  } else {
    PetscInt i, j, k, xs, ys, zs, xm, ym, zm, mx, my, mz, rstart, cstart;
    Mat B;

    ierr = DAGetInfo(daf, 0, &mz, &my, &mx, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
    ierr = DAGetCorners(daf, &zs, &ys, &xs, &zm, &ym, &xm);CHKERRQ(ierr);
    if (zs != 0) SETERRQ(1, "unexpected");
    ierr = MatCreate(((PetscObject)daf)->comm, &B);CHKERRQ(ierr);
    ierr = MatSetSizes(B, xm*ym*zm, xm*ym, mx*my*mz, mx*my);CHKERRQ(ierr);

    ierr = MatSetType(B, MATAIJ);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(B, 1, NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(B, 1, NULL, 0, NULL);CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(B, &rstart, NULL);CHKERRQ(ierr);
    ierr = MatGetOwnershipRangeColumn(B, &cstart, NULL);CHKERRQ(ierr);
    for (i=xs; i<xs+xm; i++) {
      for (j=ys; j<ys+ym; j++) {
        for (k=zs; k<zs+zm; k++) {
          PetscInt i2 = i*ym+j, i3 = i2*zm+k;
          PetscScalar val = ((k == 0 || k == mz-1) ? 0.5 : 1.) / (mz-1.); /* Integration using trapezoid rule */
          ierr = MatSetValue(B, cstart+i3, rstart+i2, val, INSERT_VALUES);CHKERRQ(ierr);
        }
      }
    }
    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatCreateMAIJ(B, sizeof(Node)/sizeof(PetscScalar), A);CHKERRQ(ierr);
    ierr = MatDestroy(B);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAGetMatrix_THI_Tridiagonal"
PetscErrorCode DAGetMatrix_THI_Tridiagonal(DA da, const MatType mtype, Mat *J)
{
  PetscErrorCode ierr;
  Mat A;
  PetscInt xm, ym, zm, dim, dof = 2, starts[3], dims[3];
  ISLocalToGlobalMapping ltog, ltogb;

  PetscFunctionBegin;
  ierr = DAGetInfo(da, &dim, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);CHKERRQ(ierr);
  if (dim != 3) SETERRQ(PETSC_ERR_ARG_WRONG, "Expected DA to be 3D");
  ierr = DAGetCorners(da, 0, 0, 0, &zm, &ym, &xm);CHKERRQ(ierr);
  ierr = DAGetISLocalToGlobalMapping(da, &ltog);CHKERRQ(ierr);
  ierr = DAGetISLocalToGlobalMappingBlck(da, &ltogb);CHKERRQ(ierr);
  ierr = MatCreate(((PetscObject)da)->comm, &A);CHKERRQ(ierr);
  ierr = MatSetSizes(A, dof*xm*ym*zm, dof*xm*ym*zm, PETSC_DETERMINE, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetType(A, mtype);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A, 6, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A, 6, PETSC_NULL, 0, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSeqBAIJSetPreallocation(A, dof, 3, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIBAIJSetPreallocation(A, dof, 3, PETSC_NULL, 0, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSeqSBAIJSetPreallocation(A, dof, 2, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPISBAIJSetPreallocation(A, dof, 2, PETSC_NULL, 0, PETSC_NULL);CHKERRQ(ierr);
  ierr = MatSetBlockSize(A, dof);CHKERRQ(ierr);
  ierr = MatSetLocalToGlobalMapping(A, ltog);CHKERRQ(ierr);
  ierr = MatSetLocalToGlobalMappingBlock(A, ltogb);CHKERRQ(ierr);
  ierr = DAGetGhostCorners(da, &starts[0], &starts[1], &starts[2], &dims[0], &dims[1], &dims[2]);CHKERRQ(ierr);
  ierr = MatSetStencil(A, dim, dims, starts, dof);CHKERRQ(ierr);
  *J = A;
  PetscFunctionReturn(0);
}

//! \brief The messy details of the THI initialization.
/*!
 * Here
 * - Lx and Ly are the domain dimensions, not the half-lengths.
 * - pism_da2 is the DA used by PISM. The 
 */
PetscErrorCode THISetup(MPI_Comm comm, DA pism_da2, PetscReal Lx, PetscReal Ly,
                        THI thi, DMMG **dmmg_out) {
  PetscErrorCode ierr;
  PetscInt i;

  // FIXME: we need to get rid of thi->alpha
  // thi->alpha = 0;

  thi->Lx = Lx;
  thi->Ly = Ly;

  ierr = DMMGCreate(PETSC_COMM_WORLD, thi->nlevels, thi, dmmg_out);CHKERRQ(ierr);
  DMMG *dmmg = *dmmg_out;;

  {
    DA da;
    int P = 2;
    ierr = PetscOptionsBegin(comm, NULL, "Grid resolution options", "");CHKERRQ(ierr);
    {
      ierr = PetscOptionsInt("-blatter_P",
                             "Number of elements in z-direction on coarse level", "", P, &P, NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd();CHKERRQ(ierr);

    // Get the parametes of the DA used in PISM so that we can create a compatible one.
    PetscInt Mx, My, mx, my;
    const PetscInt *lx, *ly;
    ierr =  DAGetInfo(pism_da2,
                      NULL,           // number of dimensions
                      &My, &Mx, NULL, // number of grid points in y, x, z directions
                      &my, &mx, NULL, // number of processors in y, x, z directions
                      NULL,           // dof
                      NULL,           // stencil width
                      NULL,           // wrap
                      NULL);          // stencil type
    CHKERRQ(ierr);
    ierr = DAGetOwnershipRanges(pism_da2, &ly, &lx, NULL); CHKERRQ(ierr);

    ierr = DACreate3d(comm, DA_YZPERIODIC, DA_STENCIL_BOX,
                      P, My, Mx, // number of grid points
                      1, my, mx, // number of processors
                      sizeof(Node)/sizeof(PetscScalar), // dof
                      1,                                // stencil width
                      &P, ly, lx,
                      &da);CHKERRQ(ierr);

    ierr = DASetFieldName(da, 0, "x-velocity");CHKERRQ(ierr);
    ierr = DASetFieldName(da, 1, "y-velocity");CHKERRQ(ierr);
    ierr = DMMGSetDM(dmmg, (DM)da);CHKERRQ(ierr);
    ierr = DADestroy(da);CHKERRQ(ierr);
  }

  if (thi->tridiagonal) {
    (DMMGGetDA(dmmg))->ops->getmatrix = DAGetMatrix_THI_Tridiagonal;
  }

  {
    /* Use the user-defined matrix type on all but the coarse level */
    ierr = DMMGSetMatType(dmmg, thi->mattype);CHKERRQ(ierr);
    /* PCREDUNDANT only works with AIJ, and so do the third-party direct
     * solvers. So when running in parallel, we can't use the faster (S)BAIJ
     * formats on the coarse level. */
    ierr = PetscFree(dmmg[0]->mtype);CHKERRQ(ierr);
    ierr = PetscStrallocpy(MATAIJ, &dmmg[0]->mtype);CHKERRQ(ierr);
  }

  /* Spectacularly ugly API, our function evaluation provides ghost values */
  ierr = PetscOptionsSetValue("-dmmg_form_function_ghost", "1");CHKERRQ(ierr);

  ierr = DMMGSetSNESLocal(dmmg, THIFunctionLocal, THIJacobianLocal_3D_Full, 0, 0);CHKERRQ(ierr);

  if (thi->tridiagonal) {
    ierr = DASetLocalJacobian(DMMGGetDA(dmmg),
                              (DALocalFunction1)THIJacobianLocal_3D_Tridiagonal);CHKERRQ(ierr);
  }

  for (i=0; i<DMMGGetLevels(dmmg); i++) {
    /* This option is only valid for the SBAIJ format. The matrices we assemble
     * are symmetric, but the SBAIJ assembly functions will complain if we
     * provide lower-triangular entries without setting this option. */
    Mat B = dmmg[i]->B;
    PetscTruth flg1, flg2;
    ierr = PetscTypeCompare((PetscObject)B, MATSEQSBAIJ, &flg1);CHKERRQ(ierr);
    ierr = PetscTypeCompare((PetscObject)B, MATMPISBAIJ, &flg2);CHKERRQ(ierr);
    if (flg1 || flg2) {
      ierr = MatSetOption(B, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  ierr = MatSetOptionsPrefix(DMMGGetB(dmmg), "thi_");CHKERRQ(ierr);
  ierr = DMMGSetFromOptions(dmmg);CHKERRQ(ierr);
  ierr = THISetDMMG(thi, dmmg);CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "DAGetPrmNodeArray"
PetscErrorCode DAGetPrmNodeArray(DA da, PrmNode ***prm)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm", (PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm_Vec", (PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm_Vec composed with given DA");
  ierr = DAVecGetArray(da2prm, X, prm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DARestorePrmNodeArray"
PetscErrorCode DARestorePrmNodeArray(DA da, PrmNode ***prm)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm", (PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm_Vec", (PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm_Vec composed with given DA");
  ierr = DAVecRestoreArray(da2prm, X, prm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DAPrmNodeArrayCommBegin(DA da)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm", (PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm_Vec", (PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm_Vec composed with given DA");
  ierr = DALocalToLocalBegin(da2prm, X, INSERT_VALUES, X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DAPrmNodeArrayCommEnd(DA da)
{
  PetscErrorCode ierr;
  DA             da2prm;
  Vec            X;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm", (PetscObject*)&da2prm);CHKERRQ(ierr);
  if (!da2prm) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm composed with given DA");
  ierr = PetscObjectQuery((PetscObject)da, "DA2Prm_Vec", (PetscObject*)&X);CHKERRQ(ierr);
  if (!X) SETERRQ(PETSC_ERR_ARG_WRONG, "No DA2Prm_Vec composed with given DA");
  ierr = DALocalToLocalEnd(da2prm, X, INSERT_VALUES, X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
