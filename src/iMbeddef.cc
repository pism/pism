// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petscda.h>
#include <cmath>
#include "iceModel.hh"
#include "extrasGSL.hh"

#if (WITH_FFTW)
#include <fftw3.h>
#endif


PetscErrorCode IceModel::createScatterToProcZero(Vec& samplep0) {
  PetscErrorCode ierr;

  ierr = DACreateNaturalVector(grid.da2, &g2natural); CHKERRQ(ierr);
  ierr = VecScatterCreateToZero(g2natural, &top0ctx, &samplep0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::destroyScatterToProcZero() {
  PetscErrorCode ierr;

  ierr = VecDestroy(g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(top0ctx); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::putOnProcZero(Vec& vlocal, Vec& onp0) {
  PetscErrorCode ierr;

  // scatter local Vec to proc zero: from a global Vec in the global ordering to 
  //    a global Vec in the natural ordering and then to a Vec on proc zero (i.e. empty on other procs)
  // requires g2, g2natural, and top0ctx to all be set up properly
  ierr = DALocalToGlobal(grid.da2,vlocal,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DAGlobalToNaturalBegin(grid.da2,g2,INSERT_VALUES,g2natural); CHKERRQ(ierr);
  ierr = DAGlobalToNaturalEnd(grid.da2,g2,INSERT_VALUES,g2natural); CHKERRQ(ierr);
  ierr = VecScatterBegin(top0ctx, g2natural,onp0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecScatterBegin(onp0,g2natural,INSERT_VALUES,SCATTER_FORWARD,top0ctx); CHKERRQ(ierr);
  ierr = VecScatterEnd(top0ctx, g2natural,onp0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecScatterEnd(onp0,g2natural,INSERT_VALUES,SCATTER_FORWARD,top0ctx); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getFromProcZero(Vec& onp0, Vec& vlocal) {
  PetscErrorCode ierr;

  // undo scatter to proc zero: put onp0 back into vlocal and communicate ghosted values
  ierr = VecScatterBegin(top0ctx, onp0,g2natural,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(top0ctx, onp0,g2natural,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = DANaturalToGlobalBegin(grid.da2,g2natural,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DANaturalToGlobalEnd(grid.da2,g2natural,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(grid.da2,g2,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(grid.da2,g2,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2,vlocal,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2,vlocal,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::conv2_same(Vec vA, const PetscInt mA, const PetscInt nA, 
                                    Vec vB, const PetscInt mB, const PetscInt nB,
                                    Vec &vresult) {
  // naively (sans FFT) convolves two *sequential* Vecs (A is mA x nA; 
  //   B is mB x nB) and returns a result Vec which is the same size as A
  // this operation is O(mA^2 nA^2)
  PetscErrorCode  ierr;
  PetscScalar     **A, **B, **result;

  ierr = VecGetArray2d(vA, mA, nA, 0, 0, &A); CHKERRQ(ierr);
  ierr = VecGetArray2d(vB, mB, nB, 0, 0, &B); CHKERRQ(ierr);
  ierr = VecGetArray2d(vresult, mA, nA, 0, 0, &result); CHKERRQ(ierr);
  for (PetscInt i=0; i < mA; i++) {
    for (PetscInt j=0; j < nA; j++) {
      PetscScalar sum = 0.0;
      for (PetscInt r = PetscMax(0, i - mB + 1); r < PetscMin(i, mA); r++) {
        for (PetscInt s = PetscMax(0, j - nB + 1); s < PetscMin(j, nA); s++) {
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


PetscErrorCode IceModel::bedDefSetup() {
  PetscErrorCode  ierr;
  const PetscInt  Mx = grid.p->Mx;
  const PetscInt  N = 2*(Mx-1);
  
  if ((doBedDef == PETSC_TRUE) && (doBedIso == PETSC_FALSE)) {
    ierr = verbPrintf(3, grid.com,
        "setting up bed deformation variables for Lingle & Clark model ...\n"); 
        CHKERRQ(ierr);
#if (WITH_FFTW==0)
    ierr = PetscPrintf(grid.com,
        "  WARNING: compiled without FFTW.  -bed_def_lc (-bed_def) will not work.\n"); 
        CHKERRQ(ierr);
#endif
    if (grid.p->Mx != grid.p->My) {  SETERRQ(1,"Mx=My required for -bed_def");  }
    if (grid.p->dx != grid.p->dy) {  SETERRQ(1,"dx=dy required for -bed_def");  }    
    
    // create scatter context (to proc zero), etc.
    ierr = createScatterToProcZero(Hstartp0); CHKERRQ(ierr);
    ierr = VecDuplicate(Hstartp0,&bedstartp0);

    // store initial H values on proc zero; note load is rho_i g (H-Hstart)
    ierr = putOnProcZero(vH,Hstartp0); CHKERRQ(ierr);
    // store initial bed values on proc zero; note bed = bedstart + plate
    ierr = putOnProcZero(vbed,bedstartp0); CHKERRQ(ierr);

    if (grid.rank == 0) {
      // initialize plate displacement
      ierr = VecCreateSeq(PETSC_COMM_SELF, N*N, &platep0fat); CHKERRQ(ierr);
      // FFT-side coefficient fields (i.e. multiplication form of operators)
      ierr = VecDuplicate(platep0fat,&vleft);
      ierr = VecDuplicate(platep0fat,&vright);

      // setup fftw stuff
#if (WITH_FFTW)
      bdin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);
      bdout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);
         
      // fftw manipulates the data in setting up a plan, so fill with nonconstant junk
      for (PetscInt i=0; i < N; i++) {
        for (PetscInt j=0; j < N; j++) {
           bdin[j + N * i][0] = static_cast<double>(i-3);
           bdin[j + N * i][1] = static_cast<double>(j*j+2);
        }
      }
      bdplanfor = fftw_plan_dft_2d(N, N, bdin, bdout, FFTW_FORWARD, FFTW_MEASURE);
      bdplanback = fftw_plan_dft_2d(N, N, bdin, bdout, FFTW_BACKWARD, FFTW_MEASURE);
#endif 
  
      // *always* initialize plate from uplift, but make sure uplift (=dbed/dt)
      // is zero if it is not actually available from data
      ierr = bed_uplift_init_lc(); CHKERRQ(ierr);

      // compare geforconv.m
      ierr = verbPrintf(3, grid.com,
         "     computing spherical elastic load response matrix ..."); CHKERRQ(ierr);
      const PetscInt  Nge = 2*Mx-1;
      ierr = VecCreateSeq(PETSC_COMM_SELF, Nge * Nge, &lrmEp0); CHKERRQ(ierr);
      PetscScalar **II;
      ierr = VecGetArray2d(lrmEp0, Nge, Nge, 0, 0, &II); CHKERRQ(ierr);
      ge_params ge_data;
      const double dx = grid.p->dx;
      const double dy = grid.p->dy;
      ge_data.dx = dx;
      ge_data.dy = dy;
      for (PetscInt i=0; i < Nge; i++) {
        for (PetscInt j=0; j < Nge; j++) {
          ge_data.p = i - Mx + 1;
          ge_data.q = j - Mx + 1;
          II[i][j] = dblquad_cubature(ge_integrand, -dx/2, dx/2, -dy/2, dy/2, 
                                      1.0e-8, &ge_data);
        }
      }
      ierr = VecRestoreArray2d(lrmEp0, Nge, Nge, 0, 0, &II); CHKERRQ(ierr);
#if 0
      ierr = PetscPrintf(grid.com, 
             "   [here is the spherical elastic LRM lrmEp0:]\n"); CHKERRQ(ierr);
      ierr = VecView(lrmEp0,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 
#endif

    } // if (grid.rank == 0)
  
    ierr = PetscPrintf(grid.com," done\n"); CHKERRQ(ierr);
  }

  if (doBedDef == PETSC_TRUE) {
    if  (doBedIso == PETSC_TRUE) {
        ierr = verbPrintf(3, grid.com,
          "setting up bed deformation vars for pointwise isostasy ...");
           CHKERRQ(ierr);
    }
    ierr = VecDuplicate(vH,&vHlast); CHKERRQ(ierr);
    ierr = VecCopy(vH,vHlast); CHKERRQ(ierr);
    ierr = VecDuplicate(vbed,&vbedlast); CHKERRQ(ierr);
    ierr = VecCopy(vbed,vbedlast); CHKERRQ(ierr);
    if  (doBedIso == PETSC_TRUE) {
      ierr = PetscPrintf(grid.com," done\n"); CHKERRQ(ierr);
    }

    lastBedDefUpdateYear = grid.p->year;
  }

  return 0;
}


PetscErrorCode IceModel::bedDefCleanup() {
  PetscErrorCode  ierr;

  if ((doBedDef == PETSC_TRUE) && (doBedIso == PETSC_FALSE)) {
    if (grid.rank == 0) {
#if (WITH_FFTW)
       fftw_destroy_plan(bdplanfor);
       fftw_destroy_plan(bdplanback);
       fftw_free(bdin);
       fftw_free(bdout);
#endif

       ierr = VecDestroy(platep0fat); CHKERRQ(ierr);
       ierr = VecDestroy(vleft); CHKERRQ(ierr);
       ierr = VecDestroy(vright); CHKERRQ(ierr);
       ierr = VecDestroy(lrmEp0); CHKERRQ(ierr);
    }  

    ierr = VecDestroy(Hstartp0); CHKERRQ(ierr);
    ierr = VecDestroy(bedstartp0); CHKERRQ(ierr);
    ierr = destroyScatterToProcZero(); CHKERRQ(ierr);
  }

  if (doBedDef == PETSC_TRUE) {
    ierr = VecDestroy(vHlast); CHKERRQ(ierr);
    ierr = VecDestroy(vbedlast); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::bedDefStepIfNeeded() {
  PetscErrorCode  ierr;

  // This is a front end to the bed deformation update system.  It updates
  // no more often than bedDefIntervalYears.

  if (doBedDef == PETSC_TRUE) {
    //ierr = PetscPrintf(grid.com, " [lastBedDefUpdateYear = %10.3f]\n",
    //           lastBedDefUpdateYear); CHKERRQ(ierr);
    // If the bed elevations are not expired, exit cleanly.
    const PetscScalar dtBedDefYears = grid.p->year - lastBedDefUpdateYear;
    if (dtBedDefYears >= bedDefIntervalYears) {
      const PetscScalar dtBedDef = dtBedDefYears * secpera;
      if (doBedIso == PETSC_TRUE) {
        ierr = bed_def_step_iso( dtBedDef ); CHKERRQ(ierr);
      } else {
        ierr = bed_def_step_lc( dtBedDef ); CHKERRQ(ierr);
      }
      // uplift = (bed - bedlast) / dt
      ierr = VecWAXPY(vuplift,-1.0,vbedlast,vbed); CHKERRQ(ierr);  
      ierr = VecScale(vuplift, 1.0 / dtBedDef); CHKERRQ(ierr); 
      // copy current values of H, bed in prep for next step
      ierr = VecCopy(vH,vHlast); CHKERRQ(ierr);
      ierr = VecCopy(vbed,vbedlast); CHKERRQ(ierr);
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
      ierr = PetscPrintf(grid.com, "b"); CHKERRQ(ierr);
      lastBedDefUpdateYear = grid.p->year;
    } else {
      ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
    }
  } else {
    ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
  }
  return 0;
}


PetscScalar IceModel::volumeChange(Vec myHp0, Vec myHstartp0) {
  // returns volume change (in m^3) since start of run 
  if (grid.rank != 0)  {
    SETERRQ(1,"volumeChange() should only be called on processor zero");
  }

  PetscErrorCode ierr;
  PetscScalar    delvolume;
  PetscScalar     **H, **Hstart;
  const PetscInt  Mx = grid.p->Mx;  // will always equal grid.p->My
  
  ierr = VecGetArray2d(myHp0, Mx, Mx, 0, 0, &H); CHKERRQ(ierr);
  ierr = VecGetArray2d(myHstartp0, Mx, Mx, 0, 0, &Hstart); CHKERRQ(ierr);
  delvolume = 0;
  const PetscScalar   a = grid.p->dx * grid.p->dy; // area unit (m^2)
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      delvolume += a * (H[i][j]-Hstart[i][j]);
    }
  }  
  ierr = VecRestoreArray2d(myHp0, Mx, Mx, 0, 0, &H); CHKERRQ(ierr);
  ierr = VecRestoreArray2d(myHstartp0, Mx, Mx, 0, 0, &Hstart); CHKERRQ(ierr);
  return delvolume;
}


PetscErrorCode IceModel::bed_def_step_iso(PetscScalar dtBedDef) {
  PetscErrorCode ierr;
  Vec vHdiff = vWork2d[0];

  const PetscScalar  f = ice.rho / bedrock.rho;
  ierr = VecCopy(vH,vHdiff); CHKERRQ(ierr);
  ierr = VecWAXPY(vHdiff,-1.0,vHlast,vH); CHKERRQ(ierr);  // Hdiff = H - Hlast
  ierr = VecWAXPY(vbed,-f,vHdiff,vbedlast); CHKERRQ(ierr);  // bed = bedlast - f (Hdiff)
  return 0;
}


PetscErrorCode IceModel::bed_uplift_init_lc() {
  // to initialize we solve: 
  //   rho_r g plate + D grad^4 plate = 0 - 2 eta |grad| uplift
  // (see Bueler, Lingle, Kallen-Brown (2006) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", sub. to Ann. Glaciol.)
  PetscErrorCode ierr;
  Vec upliftp0;

  ierr = VecDuplicate(Hstartp0, &upliftp0); CHKERRQ(ierr);
  ierr = putOnProcZero(vuplift,upliftp0); CHKERRQ(ierr);

  if (grid.rank == 0) {
#if (WITH_FFTW)
    
    PetscScalar **uplift, **plate, **left, **right;
    PetscScalar *cx;
    
    const PetscScalar pi = 3.14159265358979;
    const PetscInt    Mx = grid.p->Mx;  // will always equal grid.p->My
    const PetscInt    N = 2*(Mx-1);
    const PetscScalar L = (Mx-1) * grid.p->dx;  // L=2Lx; i.e. Z of 2; 
                                     // half-length of fat computational domain

    // spectral/FFT quantities are on fat computational grid but uplift is on thin
    ierr = VecGetArray2d(upliftp0, Mx, Mx, 0, 0, &uplift); CHKERRQ(ierr);
    ierr = VecGetArray2d(platep0fat, N, N, 0, 0, &plate); CHKERRQ(ierr);
    ierr = VecGetArray2d(vleft, N, N, 0, 0, &left); CHKERRQ(ierr);
    ierr = VecGetArray2d(vright, N, N, 0, 0, &right); CHKERRQ(ierr);

    // fft2(uplift)
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        bdin[j + N * i][0] = 0.0;
        bdin[j + N * i][1] = 0.0;
      }
    }
    for (PetscInt i=0; i < Mx; i++) {
      for (PetscInt j=0; j < Mx; j++) {
        bdin[( j+(Mx-1)/2 ) + N * ( i+(Mx-1)/2 )][0] = uplift[i][j];
      }
    }
    fftw_execute(bdplanfor);

    // Matlab version:  cx=(pi/L)*[0:N/2 N/2-1:-1:1]
    cx = new PetscScalar[N];
    for (PetscInt i=0; i <= N/2; i++) {
      cx[i] = (pi/L) * (PetscScalar) i;
    }
    for (PetscInt i=N/2+1; i < N; i++) {
      cx[i] = (pi/L) * (PetscScalar) (N - i);
    }
    // compute left and right coefficients
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        const PetscScalar cclap = cx[i]*cx[i] + cx[j]*cx[j];
        left[i][j] = bedrock.rho * grav + bedrock.D * cclap * cclap;
        right[i][j] = -2 * bedrock.eta * sqrt(cclap);
      }
    }
    delete [] cx;

    // Matlab version would look something like:
    //        frhs = right.*fft2(uplift);
    //        u = real(ifft2( frhs./left ));
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        bdin[j + N * i][0] = (right[i][j] * bdout[j + N * i][0]) / left[i][j];
        bdin[j + N * i][1] = (right[i][j] * bdout[j + N * i][1]) / left[i][j];
      }
    }
    fftw_execute(bdplanback);
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        plate[i][j] = bdout[j + N * i][0] / ((PetscScalar) (N*N));
      }
    }
    ierr = VecRestoreArray2d(vleft, N, N, 0, 0, &left); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(vright, N, N, 0, 0, &right); CHKERRQ(ierr);
    
    ierr = VecRestoreArray2d(platep0fat, N, N, 0, 0, &plate); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(upliftp0, Mx, Mx, 0, 0, &uplift); CHKERRQ(ierr);
#endif
  }

  ierr = VecDestroy(upliftp0); CHKERRQ(ierr);  
  return 0;
}


PetscErrorCode IceModel::bed_def_step_lc(PetscScalar dtBedDef) {
  PetscErrorCode ierr;
  Vec Hp0, Hdiffp0, bedp0;

  ierr = VecDuplicate(Hstartp0, &bedp0); CHKERRQ(ierr); // merely create space

  ierr = VecDuplicate(Hstartp0, &Hp0); CHKERRQ(ierr);
  ierr = putOnProcZero(vH,Hp0); CHKERRQ(ierr);
  ierr = VecDuplicate(Hstartp0, &Hdiffp0); CHKERRQ(ierr);
  ierr = VecWAXPY(Hdiffp0,-1.0,Hstartp0,Hp0); CHKERRQ(ierr);  // Hdiff = H - Hstart

  if (grid.rank == 0) {
#if (WITH_FFTW)
    // see Bueler, Lingle, Kallen-Brown (2006) "Fast computation of a 
    // viscoelastic deformable Earth model for ice sheet simulations", sub. to 
    // Ann. Glaciol.
    
    PetscScalar **Hdiff, **bed, **bedstart, **plate, **left, **right;
    PetscScalar *cx;
    const PetscScalar pi = 3.14159265358979;
    const PetscInt    Mx = grid.p->Mx;  // will always equal grid.p->My
    const PetscInt    N = 2*(Mx-1);
    const PetscScalar L = (Mx-1) * grid.p->dx;  // L=2Lx; i.e. Z of 2; 
                                     // half-length of fat computational domain

    // note ice thicknesses and bed elevations only on physical ("thin") grid
    //   while spectral/FFT quantities are on fat computational grid
    ierr = VecGetArray2d(platep0fat, N, N, 0, 0, &plate); CHKERRQ(ierr);
    ierr = VecGetArray2d(vleft, N, N, 0, 0, &left); CHKERRQ(ierr);
    ierr = VecGetArray2d(vright, N, N, 0, 0, &right); CHKERRQ(ierr);

    // Matlab: sszz=-rhoi*g*H;  (here H -> H-Hstart) 
    //         COMPUTE fft2(dt*sszz);
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        bdin[j + N * i][0] = 0.0;
        bdin[j + N * i][1] = 0.0;
      }
    }
    ierr = VecGetArray2d(Hdiffp0, Mx, Mx, 0, 0, &Hdiff); CHKERRQ(ierr);
    for (PetscInt i=0; i < Mx; i++) {
      for (PetscInt j=0; j < Mx; j++) {
        const PetscScalar sszz = - ice.rho * grav * Hdiff[i][j];
        bdin[( j+(Mx-1)/2 ) + N * ( i+(Mx-1)/2 )][0] = dtBedDef * sszz;
      }
    }
    ierr = VecRestoreArray2d(Hdiffp0, Mx, Mx, 0, 0, &Hdiff); CHKERRQ(ierr);
    fftw_execute(bdplanfor);
    fftw_complex  *loadhat;  // actually a 2D array but must be flat for fftw
    loadhat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N);
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        loadhat[j + N * i][0] = bdout[j + N * i][0];
        loadhat[j + N * i][1] = bdout[j + N * i][1];
      }
    }

    //         COMPUTE fft2(uun);
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        bdin[j + N * i][0] = plate[i][j];
        bdin[j + N * i][1] = 0.0;
      }
    }
    fftw_execute(bdplanfor);

    // Matlab version:  cx=(pi/L)*[0:N/2 N/2-1:-1:1];
    cx = new PetscScalar[N];
    for (PetscInt i=0; i <= N/2; i++) {
      cx[i] = (pi/L) * (PetscScalar) i;
    }
    for (PetscInt i=N/2+1; i < N; i++) {
      cx[i] = (pi/L) * (PetscScalar) (N - i);
    }
    // compute left and right coefficients; note they depend on actual time step
    //   and thus cannot be precomputed
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        const PetscScalar cclap = cx[i]*cx[i] + cx[j]*cx[j];
        const PetscScalar part1 = 2 * bedrock.eta * sqrt(cclap);
        const PetscScalar part2 = (dtBedDef/2) *
                      (bedrock.rho * grav + bedrock.D * cclap * cclap);
        left[i][j] = part1 + part2;
        right[i][j] = part1 - part2;
      }
    }
    delete [] cx;

    //         frhs=right.*fft2(uun) + fft2(dt*sszz);
    //         uun1=real(ifft2( frhs./left ));
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        bdin[j + N * i][0] = (right[i][j] * bdout[j + N * i][0] + 
                                 loadhat[j + N * i][0]) / left[i][j];
        bdin[j + N * i][1] = (right[i][j] * bdout[j + N * i][1] + 
                                 loadhat[j + N * i][1]) / left[i][j];
      }
    }
    fftw_execute(bdplanback);
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        plate[i][j] = bdout[j + N * i][0] / ((PetscScalar) (N*N));
      }
    }
    fftw_free(loadhat);
    ierr = VecRestoreArray2d(vleft, N, N, 0, 0, &left); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(vright, N, N, 0, 0, &right); CHKERRQ(ierr);

    // tweak: find average value along "distant" bdry of [-ZL,ZL]X[-ZL,ZL]
    // (will remove it:   uun1=uun1-( sum(uun1(1,:))+sum(uun1(:,1)) )/(2*N);)
    PetscScalar   sum = 0.0;
    for (PetscInt i=0; i < N; i++) {
      sum += plate[i][0] + plate[0][i];
    }
    // tweak cont: replace far field with value for equiv disc load which has R0=Lx*(2/3)=L/3 
    // (instead of 1000km in Matlab code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2);  % trapezoid rule)
    const PetscScalar Requiv = L/3.0;
    const PetscScalar Hequiv = volumeChange(Hp0,Hstartp0) / (pi * Requiv * Requiv);
    // note: viscDisc() returns 0.0 if GSL not available so effect is to just remove far field value
    // FIXME: is this the right time if model is restarted from .pb file?
    const PetscScalar vd_time = (grid.p->year - startYear) * secpera;
    const PetscScalar discshift = viscDisc(vd_time,Hequiv,Requiv,L,      
                                           bedrock.rho,grav,bedrock.D,bedrock.eta)
                                  - sum/(2*N); 
    for (PetscInt i=0; i < N; i++) {
      for (PetscInt j=0; j < N; j++) {
        plate[i][j] += discshift;
      }
    }
#if 0
    ierr = PetscPrintf(grid.com, "   [Hequiv = %12.4f]\n",Hequiv); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, "   [boundary average = %12.4f]\n",sum/(2*N)); CHKERRQ(ierr);
    ierr = PetscPrintf(grid.com, "   [viscDisc() result = %12.4f]\n",discshift+sum/(2*N)); CHKERRQ(ierr);
#endif

    // Matlab:     ue=rhoi*conv2(H,II,'same');  (but H --> H-Hstart here)
    // first, bed = ue;
#if 1
    const PetscInt  Nge = 2*Mx-1;
    ierr = conv2_same(Hdiffp0,Mx,Mx,lrmEp0,Nge,Nge,bedp0);  CHKERRQ(ierr);
#else
    ierr = VecSet(bedp0,0.0); CHKERRQ(ierr);
#endif
    ierr = VecScale(bedp0,ice.rho);  CHKERRQ(ierr);

    // (new bed) = ue + (bed start) + plate; uses only central part of plate
    ierr = VecGetArray2d(bedp0, Mx, Mx, 0, 0, &bed); CHKERRQ(ierr);
    ierr = VecGetArray2d(bedstartp0, Mx, Mx, 0, 0, &bedstart); CHKERRQ(ierr);    
    for (PetscInt i=0; i < Mx; i++) {
      for (PetscInt j=0; j < Mx; j++) {
        bed[i][j] += bedstart[i][j] + plate[i+(Mx-1)/2][j+(Mx-1)/2];
      }
    }
    ierr = VecRestoreArray2d(bedp0, Mx, Mx, 0, 0, &bed); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(bedstartp0, Mx, Mx, 0, 0, &bedstart); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(platep0fat, N, N, 0, 0, &plate); CHKERRQ(ierr);
#endif
  }

#if 0
  ierr = PetscPrintf(grid.com, "   [here is bedp0 AGAIN:]\n"); CHKERRQ(ierr);
  ierr = VecView(bedp0,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); 
#endif

  // undo scatter of bed at beginning; put bedp0 back into vbed;
  // communicates ghosted values of bed (needed because they effect surface slope)
  // (note Hp0 is discarded, which is fine)
  ierr = getFromProcZero(bedp0,vbed); CHKERRQ(ierr);

  ierr = VecDestroy(Hp0); CHKERRQ(ierr);
  ierr = VecDestroy(Hdiffp0); CHKERRQ(ierr);
  ierr = VecDestroy(bedp0); CHKERRQ(ierr);
  
  return 0;
}
