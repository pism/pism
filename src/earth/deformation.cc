// Copyright (C) 2004-2009, 2011 Ed Bueler
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
#include <petscvec.h>
#include "pism_const.hh"
#include "matlablike.hh"
#include "greens.hh"

#if (PISM_HAVE_FFTW)
#include <fftw3.h>
#endif

#include "deformation.hh"


BedDeformLC::BedDeformLC() {
  settingsDone = PETSC_FALSE;
  allocDone = PETSC_FALSE;
}


BedDeformLC::~BedDeformLC() {
  if (allocDone == PETSC_TRUE) {
#if (PISM_HAVE_FFTW)
    fftw_destroy_plan(bdplanfor);
    fftw_destroy_plan(bdplanback);
    fftw_free(bdin);
    fftw_free(bdout);
#endif
    VecDestroy(Hdiff);
    VecDestroy(dbedElastic);
    VecDestroy(platefat);
    VecDestroy(plateoffset);
    VecDestroy(vleft);
    VecDestroy(vright);
    VecDestroy(lrmE);
    delete [] cx;  delete [] cy;
  }
  allocDone = PETSC_FALSE;
}


PetscErrorCode BedDeformLC::settings(const NCConfigVariable &config,
                  PetscTruth  myinclude_elastic,
                  PetscInt myMx, PetscInt myMy, PetscScalar mydx, PetscScalar mydy,
                  PetscInt myZ, PetscScalar myicerho,
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
  rho    = config.get("lithosphere_density");
  eta    = config.get("mantle_viscosity");
  D      = config.get("lithosphere_flexural_rigidity");

  standard_gravity = config.get("standard_gravity");

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
      lft[i][j] = rho * standard_gravity + D * cclap * cclap;
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
      const PetscScalar sszz = - icerho * standard_gravity * dH[i][j];
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
      const PetscScalar part2 = (dt/2.0) * (rho * standard_gravity + D * cclap * cclap);
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
  const PetscScalar discshift = viscDisc(vd_time,Hequiv,Requiv,Lav,rho,standard_gravity,D,eta) - av; 
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

