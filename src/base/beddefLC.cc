// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#if (PISM_HAVE_FFTW)
#include <fftw3.h>
#endif
#include <petscvec.h>
#include "pism_const.hh"
#include "../num/extrasGSL.hh"
#include "beddefLC.hh"


double ge_integrand(unsigned ndimMUSTBETWO, const double* xiANDeta, void* paramsIN) {
  // Matlab:  function z=integrand(xi,eta,dx,dy,p,q)

  if (ndimMUSTBETWO != 2) { perror("ge_integrand only defined for 2 variables"); }
  
  // data here is from Lingle & Clark (1985)
  double rmkm[42] =
    {0.0, 0.011,  0.111,  1.112,  2.224,  3.336,  4.448,  6.672,  8.896,  11.12, 17.79,
          22.24,  27.80,  33.36,  44.48,  55.60,  66.72,  88.96,  111.2,  133.4, 177.9,
          222.4,  278.0,  333.6,  444.8,  556.0,  667.2,  778.4,  889.6, 1001.0, 1112.0,
         1334.0, 1779.0, 2224.0, 2780.0, 3336.0, 4448.0, 5560.0, 6672.0, 7784.0, 8896.0,
        10008.0};
  // rm = rmkm * 1e3 (remember to convert to meters); GE /(10^12 rm) is
  //    vertical displacement in meters
  // (GE(r=0) has been computed by linear extrapolation:  GE(0) := -33.6488)
  double GE[42] =
    {-33.6488, -33.64, -33.56, -32.75, -31.86, -30.98, -30.12, -28.44, -26.87, -25.41,
               -21.80, -20.02, -18.36, -17.18, -15.71, -14.91, -14.41, -13.69, -13.01,
               -12.31, -10.95, -9.757, -8.519, -7.533, -6.131, -5.237, -4.660, -4.272,
               -3.999, -3.798, -3.640, -3.392, -2.999, -2.619, -2.103, -1.530, -0.292,
                0.848,  1.676,  2.083,  2.057,  1.643};  

  struct ge_params*  params = (struct ge_params*) paramsIN;
  const double  dx = params->dx;
  const double  dy = params->dy;
  const int     p = params->p;
  const int     q = params->q;
  const double  xishift = (double) p * dx - xiANDeta[0];
  const double  etashift = (double) q * dy - xiANDeta[1];
  const double  r = sqrt(xishift * xishift + etashift * etashift);
  
  double  z;
  if (r < 0.01) {
    z = GE[0]/(rmkm[1] * 1.0e3 * 1.0e12);
  } else if (r > rmkm[41] * 1.0e3) {
    z = 0.0;
  } else {
    z = interp1_linear(rmkm,GE,42,r / 1.0e3) / (r * 1.0e12);
  }
  return z;
}


double viscDiscIntegrand (double kap, void * paramsIN) {
// Matlab:  function y=integrand(kap,rg,D,t,eta,R0,rk)
//            beta=rg + D*kap.^4;
//            expdiff=exp(-beta*t./(2*eta*kap))-ones(size(kap));
//            y=expdiff.*besselj(1.0,kap*R0).*besselj(0.0,kap*rk)./beta;

  struct vd_params*  params = (struct vd_params*) paramsIN;
  const double       t = params->t;
  const double       R0 = params->R0;
  const double       rk = params->rk; 
  const double       rho = params->rho; 
  const double       grav = params->grav; 
  const double       D = params->D;
  const double       eta = params->eta;
  const double       beta = rho * grav + D * pow(kap,4.0);
  const double       expdiff = exp(-beta * t / (2.0 * eta * kap)) - 1.0;
  return expdiff * gsl_sf_bessel_J1(kap * R0) * gsl_sf_bessel_J0(kap * rk) / beta;
}


double viscDisc(double t, double H0, double R0, double r, 
                double rho, double grav, double D, double eta) {
  // t in seconds; H0, R0, r in meters

  const double      ABSTOL = 1.0e-10;
  const double      RELTOL = 1.0e-14;
  const int         N_gsl_workspace = 1000;
  gsl_integration_workspace*
                    w = gsl_integration_workspace_alloc(N_gsl_workspace);
  double*           pts;
  const int         lengthpts = 142;
  
  // Matlab:  pts=[10.^(-3:-0.05:-10) 1.0e-14];
  pts = new double[lengthpts];
  for (int j=0; j < lengthpts-1; j++) {
    pts[j] = pow(10.0,-3.0 - 0.05 * (double) j);
  }
  pts[lengthpts-1] = 1.0e-14;

  // result=quadl(@integrand,pts(1),100.0*pts(1),TOL,0,rg,D,t,eta,R0,rk); % kap->infty tail
  gsl_function      F;
  struct vd_params  params = { t, R0, r, rho, grav, D, eta };
  double            result, error;
  F.function = &viscDiscIntegrand;
  F.params = &params;
  // regarding tolerance: request is for convergence of all digits and relative tolerance RELTOL
  gsl_integration_qag (&F, pts[1], 100.0*pts[1], ABSTOL, RELTOL, N_gsl_workspace, 
                       GSL_INTEG_GAUSS21, w, &result, &error);

  double  sum = result; 
  // for j=1:length(pts)-1
  //   result=result+quadl(@integrand,pts(j+1),pts(j),TOL,0,rg,D,t,eta,R0,rk);
  // end
  for (int j=0; j < lengthpts-1; j++) {
    gsl_integration_qag (&F, pts[j+1], pts[j], ABSTOL, RELTOL, N_gsl_workspace, 
                         GSL_INTEG_GAUSS21, w, &result, &error);
    sum += result;
  }
  
  delete [] pts;
  gsl_integration_workspace_free(w);
  // u(k)=rhoi*g*H0*R0*result;
  return rho * grav * H0 * R0 * sum;
}


PetscErrorCode conv2_same(Vec vA, const PetscInt mA, const PetscInt nA, 
                                    Vec vB, const PetscInt mB, const PetscInt nB,
                                    Vec &vresult) {
  // naively (sans FFT) convolves two *sequential* Vecs (A is mA x nA; 
  // B is mB x nB) and returns a result Vec which is the same size as A
  // effectively pads by zero, in that the result matches the discrete but
  // infinite case case where A(i,j) and B(i,j) are defined for
  // -\infty < i,j < \infty but A(i,j)=0 if i<0 or i>mA-1 or j<0 or j>nA-1
  // and B(i,j)=0 if i<0 or i>mB-1 or j<0 or j>nB-1
  // this operation is O(mA^2 nA^2), but an alternate FFT implementation would
  // give O(m n log(m n)), presumably
  PetscErrorCode  ierr;
  PetscScalar     **A, **B, **result;

  ierr = VecGetArray2d(vA, mA, nA, 0, 0, &A); CHKERRQ(ierr);
  ierr = VecGetArray2d(vB, mB, nB, 0, 0, &B); CHKERRQ(ierr);
  ierr = VecGetArray2d(vresult, mA, nA, 0, 0, &result); CHKERRQ(ierr);
  for (PetscInt i=0; i < mA; i++) {
    for (PetscInt j=0; j < nA; j++) {
      PetscScalar sum = 0.0;
      for (PetscInt r = PetscMax(0, i - mB + 1); r < PetscMin(mA, i); r++) {
        for (PetscInt s = PetscMax(0, j - nB + 1); s < PetscMin(nA, j); s++) {
          sum += A[r][s] * B[i - r][j - s];
        }
      }
      result[i][j] = sum;
    }
  }
  ierr = VecRestoreArray2d(vA, mA, nA, 0, 0, &A); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(vB, mB, nB, 0, 0, &B); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(vresult, mA, nA, 0, 0, &result); CHKERRQ(ierr);
  return 0;
}


BedDeformLC::BedDeformLC() {
  settingsDone = PETSC_FALSE;
  allocDone = PETSC_FALSE;
}


BedDeformLC::~BedDeformLC() {
  if (allocDone == PETSC_TRUE) {
    PetscErrorCode  ierr;
#if (PISM_HAVE_FFTW)
    fftw_destroy_plan(bdplanfor);
    fftw_destroy_plan(bdplanback);
    fftw_free(bdin);
    fftw_free(bdout);
#endif
    ierr = VecDestroy(Hdiff);
    ierr = VecDestroy(dbedElastic);
    ierr = VecDestroy(platefat);
    ierr = VecDestroy(plateoffset);
    ierr = VecDestroy(vleft);
    ierr = VecDestroy(vright);
    ierr = VecDestroy(lrmE);
    delete [] cx;  delete [] cy;
  }
  allocDone = PETSC_FALSE;
}


PetscErrorCode BedDeformLC::settings(
                  PetscTruth  myinclude_elastic,
                  PetscInt myMx, PetscInt myMy, PetscScalar mydx, PetscScalar mydy,
                  PetscInt myZ, PetscScalar myicerho,
                  PetscScalar myrho, PetscScalar myeta, PetscScalar myD,
                  Vec* myHstart, Vec* mybedstart, Vec* myuplift, 
                  Vec* myH, Vec* mybed) {

  // set parameters
  include_elastic = myinclude_elastic;
  Mx     = myMx;
  My     = myMy;
  dx     = mydx;
  dy     = mydy;
  Z      = myZ; 
  icerho = myicerho;
  rho    = myrho;
  eta    = myeta;
  D      = myD;

  // derive more parameters
  Lx     = ((Mx-1)/2) * dx;
  Ly     = ((My-1)/2) * dy; 
  Nx     = Z*(Mx-1);
  Ny     = Z*(My-1);      
  fatLx  = (Nx/2) * dx;
  fatLy  = (Ny/2) * dy; 
  Nxge   = Nx + 1;
  Nyge   = Ny + 1;
  i0_plate = (Z-1)*(Mx-1)/2;
  j0_plate = (Z-1)*(My-1)/2;
  
  // attach to existing (must be allocated!) sequential Vecs
  H        = myH;
  bed      = mybed;
  Hstart   = myHstart;
  bedstart   = mybedstart;
  uplift   = myuplift;

  settingsDone = PETSC_TRUE;
  return 0;
}


PetscErrorCode BedDeformLC::alloc() {
  PetscErrorCode  ierr;
  if (settingsDone == PETSC_FALSE) {
    SETERRQ(1,"BedDeformLC must be set with settings() before alloc()\n");
  }
  if (allocDone == PETSC_TRUE) {
    SETERRQ(2,"BedDeformLC already allocated\n");
  }
#if (PISM_HAVE_FFTW==0)
  SETERRQ(3,"BedDeformLC will not work without FFTW");
#endif

  ierr = VecDuplicate(*H,&Hdiff); CHKERRQ(ierr);  // allocate working space
  ierr = VecDuplicate(*H,&dbedElastic); CHKERRQ(ierr);  // allocate working space
  
  // allocate plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, Nx * Ny, &platefat); CHKERRQ(ierr);
  ierr = VecDuplicate(platefat,&plateoffset); CHKERRQ(ierr);
  // FFT-side coefficient fields (i.e. multiplication form of operators)
  ierr = VecDuplicate(platefat,&vleft); CHKERRQ(ierr);
  ierr = VecDuplicate(platefat,&vright); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, Nxge * Nyge, &lrmE); CHKERRQ(ierr);

  // setup fftw stuff: FFTW builds "plans" based on observed performance
#if (PISM_HAVE_FFTW)
  bdin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  bdout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
         
  // fftw manipulates the data in setting up a plan, so fill with nonconstant junk
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
       bdin[j + Ny * i][0] = static_cast<double>(i-3);
       bdin[j + Ny * i][1] = static_cast<double>(j*j+2);
    }
  }
  bdplanfor  = fftw_plan_dft_2d(Nx, Ny, bdin, bdout, FFTW_FORWARD, FFTW_MEASURE);
  bdplanback = fftw_plan_dft_2d(Nx, Ny, bdin, bdout, FFTW_BACKWARD, FFTW_MEASURE);
#endif 
  
  // coeffs for Fourier spectral method Laplacian
  // Matlab version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  cx = new PetscScalar[Nx];
  cy = new PetscScalar[Ny];
  for (PetscInt i=0; i <= Nx/2; i++) {
    cx[i] = (pi/fatLx) * (PetscScalar) i;
  }
  for (PetscInt i=Nx/2+1; i < Nx; i++) {
    cx[i] = (pi/fatLx) * (PetscScalar) (Nx - i);
  }
  for (PetscInt j=0; j <= Ny/2; j++) {
    cy[j] = (pi/fatLy) * (PetscScalar) j;
  }
  for (PetscInt j=Ny/2+1; j < Ny; j++) {
    cy[j] = (pi/fatLy) * (PetscScalar) (Ny - j);
  }

  // compare geforconv.m
  if (include_elastic == PETSC_TRUE) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
           "     computing spherical elastic load response matrix ..."); CHKERRQ(ierr);
    PetscScalar **II;
    ierr = VecGetArray2d(lrmE, Nxge, Nyge, 0, 0, &II); CHKERRQ(ierr);
    ge_params ge_data;
    ge_data.dx = dx;
    ge_data.dy = dy;
    //const PetscInt imid_ge = Nx/2, jmid_ge = Ny/2;
    for (PetscInt i=0; i < Nxge; i++) {
      for (PetscInt j=0; j < Nyge; j++) {
        ge_data.p = i;
        ge_data.q = j;
        //ge_data.p = i - imid_ge;
        //ge_data.q = j - jmid_ge;
        II[i][j] = dblquad_cubature(ge_integrand, -dx/2, dx/2, -dy/2, dy/2, 
                                    1.0e-8, &ge_data);
      }
    }
    ierr = VecRestoreArray2d(lrmE, Nxge, Nyge, 0, 0, &II); CHKERRQ(ierr);
    // ierr = VecView(lrmE,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF," done\n"); CHKERRQ(ierr);
  }

  allocDone = PETSC_TRUE;
  return 0;
}


PetscErrorCode BedDeformLC::uplift_init() {
  // to initialize we solve: 
  //   rho_r g U + D grad^4 U = 0 - 2 eta |grad| uplift
  // where U=plateoffset; yes it really should have "0" on right side because
  // load for future times will be "-rho g (H-Hstart)", which is zero if no geometry 
  // change
  // compare equation (16) in 
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105
  // [NOTE PROBABLE SIGN ERROR in eqn (16)?:  load "rho g H" should be "- rho g H"]
    
#if (PISM_HAVE_FFTW)
  PetscErrorCode ierr;
  PetscScalar **upl, **pla, **lft, **rgt;

  // spectral/FFT quantities are on fat computational grid but uplift is on thin
  ierr = VecGetArray2d(*uplift, Mx, My, 0, 0, &upl); CHKERRQ(ierr);
  ierr = VecGetArray2d(plateoffset, Nx, Ny, 0, 0, &pla); CHKERRQ(ierr);
  ierr = VecGetArray2d(vleft, Nx, Ny, 0, 0, &lft); CHKERRQ(ierr);
  ierr = VecGetArray2d(vright, Nx, Ny, 0, 0, &rgt); CHKERRQ(ierr);

  // fft2(uplift)
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      bdin[j + Ny * i][0] = 0.0;
      bdin[j + Ny * i][1] = 0.0;
    }
  }
  for (PetscInt i=0; i < Mx; i++) {
    for (PetscInt j=0; j < My; j++) {
      bdin[(j + j0_plate) + Ny * (i + i0_plate)][0] = upl[i][j];
    }
  }
  fftw_execute(bdplanfor);

  // compute left and right coefficients
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      const PetscScalar cclap = cx[i]*cx[i] + cy[j]*cy[j];
      lft[i][j] = rho * grav + D * cclap * cclap;
      rgt[i][j] = -2.0 * eta * sqrt(cclap);
    }
  }

  // Matlab version:
  //        frhs = right.*fft2(uplift);
  //        u = real(ifft2( frhs./left ));
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      bdin[j + Ny * i][0] = (rgt[i][j] * bdout[j + Ny * i][0]) / lft[i][j];
      bdin[j + Ny * i][1] = (rgt[i][j] * bdout[j + Ny * i][1]) / lft[i][j];
    }
  }
  fftw_execute(bdplanback);
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      pla[i][j] = bdout[j + Ny * i][0] / ((PetscScalar) (Nx * Ny));
    }
  }

  ierr = VecRestoreArray2d(*uplift, Mx, My, 0, 0, &upl); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(plateoffset, Nx, Ny, 0, 0, &pla); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(vleft, Nx, Ny, 0, 0, &lft); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(vright, Nx, Ny, 0, 0, &rgt); CHKERRQ(ierr);

  ierr = VecCopy(plateoffset,platefat); CHKERRQ(ierr);
#endif

  return 0;
}


PetscErrorCode BedDeformLC::step(const PetscScalar dtyear, const PetscScalar yearFromStart) {
  // solves: 
  //     (2 eta |grad| U^{n+1}) + (dt/2) * ( rho_r g U^{n+1} + D grad^4 U^{n+1} )
  //   = (2 eta |grad| U^n) - (dt/2) * ( rho_r g U^n + D grad^4 U^n ) - dt * rho g Hstart
  // where U=plate; see equation (7) in 
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105
  PetscErrorCode ierr;
  PetscScalar **b, **bs, **dbE, **pla, **plaoff;

  // update Hdiff from H
  ierr = VecWAXPY(Hdiff, -1, *Hstart, *H); CHKERRQ(ierr);

#if (PISM_HAVE_FFTW)
  const PetscScalar dt = dtyear * secpera;
  PetscScalar **dH, **lft, **rgt;

  // note ice thicknesses and bed elevations only on physical ("thin") grid
  //   while spectral/FFT quantities are on fat computational grid
  ierr = VecGetArray2d(platefat, Nx, Ny, 0, 0, &pla); CHKERRQ(ierr);
  ierr = VecGetArray2d(vleft, Nx, Ny, 0, 0, &lft); CHKERRQ(ierr);
  ierr = VecGetArray2d(vright, Nx, Ny, 0, 0, &rgt); CHKERRQ(ierr);

  // Matlab: sszz=-rhoi*g*H;  (here H -> H-Hstart) 
  //         COMPUTE fft2(dt*sszz);
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      bdin[j + Ny * i][0] = 0.0;
      bdin[j + Ny * i][1] = 0.0;
    }
  }
  ierr = VecGetArray2d(Hdiff, Mx, My, 0, 0, &dH); CHKERRQ(ierr);
  for (PetscInt i=0; i < Mx; i++) {
    for (PetscInt j=0; j < My; j++) {
      const PetscScalar sszz = - icerho * grav * dH[i][j];
      bdin[(j + j0_plate) + Ny * (i + i0_plate)][0] = dt * sszz;
    }
  }
  ierr = VecRestoreArray2d(Hdiff, Mx, My, 0, 0, &dH); CHKERRQ(ierr);
  fftw_execute(bdplanfor);
  fftw_complex  *loadhat;  // actually a 2D array but must be flat for fftw
  loadhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      loadhat[j + Ny * i][0] = bdout[j + Ny * i][0];
      loadhat[j + Ny * i][1] = bdout[j + Ny * i][1];
    }
  }

  //         COMPUTE fft2(uun);
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      bdin[j + Ny * i][0] = pla[i][j];
      bdin[j + Ny * i][1] = 0.0;
    }
  }
  fftw_execute(bdplanfor);

  // compute left and right coefficients; note they depend on actual time step
  //   and thus cannot be precomputed
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      const PetscScalar cclap = cx[i]*cx[i] + cy[j]*cy[j];
      const PetscScalar part1 = 2.0 * eta * sqrt(cclap);
      const PetscScalar part2 = (dt/2.0) * (rho * grav + D * cclap * cclap);
      lft[i][j] = part1 + part2;
      rgt[i][j] = part1 - part2;
    }
  }

  //         frhs=right.*fft2(uun) + fft2(dt*sszz);
  //         uun1=real(ifft2( frhs./left ));
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      bdin[j + Ny * i][0] = (rgt[i][j] * bdout[j + Ny * i][0] + 
                               loadhat[j + Ny * i][0]) / lft[i][j];
      bdin[j + Ny * i][1] = (rgt[i][j] * bdout[j + Ny * i][1] + 
                               loadhat[j + Ny * i][1]) / lft[i][j];
    }
  }
  fftw_execute(bdplanback);
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      pla[i][j] = bdout[j + Ny * i][0] / ((PetscScalar) (Nx * Ny));
    }
  }
  fftw_free(loadhat);

  ierr = VecRestoreArray2d(vleft, Nx, Ny, 0, 0, &lft); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(vright, Nx, Ny, 0, 0, &rgt); CHKERRQ(ierr);
#endif

  // now tweak
#if 1
  // find average value along "distant" bdry of [-fatLx,fatLx]X[-fatLy,fatLy]
  // note domain is periodic, so think of cut locus of torus (!)
  // (will remove it:   uun1=uun1-( sum(uun1(1,:))+sum(uun1(:,1)) )/(2*N);)
  PetscScalar   av = 0.0;
  for (PetscInt i=0; i < Nx; i++) {
    av += pla[i][0];
  }
  for (PetscInt j=0; j < Ny; j++) {
    av += pla[0][j];
  }
  av = av / ((PetscScalar) (Nx + Ny));
  // tweak cont: replace far field with value for equiv disc load which has R0=Lx*(2/3)=L/3 
  // (instead of 1000km in Matlab code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2);  % trapezoid rule)
  const PetscScalar Lav = (fatLx + fatLy)/2.0;
  const PetscScalar Requiv = Lav * (2.0/3.0);
  PetscScalar delvolume;
  ierr = VecSum(Hdiff, &delvolume); CHKERRQ(ierr);
  delvolume = delvolume * dx * dy;  // make into a volume
  const PetscScalar Hequiv = delvolume / (pi * Requiv * Requiv);

  const PetscScalar vd_time = yearFromStart * secpera;
  const PetscScalar discshift = viscDisc(vd_time,Hequiv,Requiv,Lav,rho,grav,D,eta) - av; 
  for (PetscInt i=0; i < Nx; i++) {
    for (PetscInt j=0; j < Ny; j++) {
      pla[i][j] += discshift;
    }
  }
#endif

  // now compute elastic response if desired; bed = ue at end of this block
  if (include_elastic == PETSC_TRUE) {
    // Matlab:     ue=rhoi*conv2(H-Hstart,II,'same')
    ierr = conv2_same(Hdiff,Mx,My,lrmE,Nxge,Nyge,dbedElastic);  CHKERRQ(ierr);
    ierr = VecScale(dbedElastic,icerho);  CHKERRQ(ierr);
  } else {
    ierr = VecSet(dbedElastic,0.0); CHKERRQ(ierr);
  }

  // now sum contributions to get new bed elevation:
  //    (new bed) = ue + (bed start) + plate
  // (but use only central part of plate if Z>1)
  ierr = VecGetArray2d(plateoffset, Nx, Ny, 0, 0, &plaoff); CHKERRQ(ierr);
  ierr = VecGetArray2d(*bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
  ierr = VecGetArray2d(*bedstart, Mx, My, 0, 0, &bs); CHKERRQ(ierr);
  ierr = VecGetArray2d(dbedElastic, Mx, My, 0, 0, &dbE); CHKERRQ(ierr);
  for (PetscInt i=0; i < Mx; i++) {
    for (PetscInt j=0; j < My; j++) {
      b[i][j] = bs[i][j] + dbE[i][j] // add elastic first
                 + (pla[i + i0_plate][j + j0_plate] - plaoff[i + i0_plate][j + j0_plate]);
    }
  }
  ierr = VecRestoreArray2d(*bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(*bedstart, Mx, My, 0, 0, &bs); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(dbedElastic, Mx, My, 0, 0, &dbE); CHKERRQ(ierr);    
  ierr = VecRestoreArray2d(plateoffset, Nx, Ny, 0, 0, &plaoff); CHKERRQ(ierr);

  ierr = VecRestoreArray2d(platefat, Nx, Ny, 0, 0, &pla); CHKERRQ(ierr);
  return 0;
}

