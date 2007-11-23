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
#include <petscda.h>
#include "../exact/exactTestsFG.h" 
#include "../exact/exactTestK.h" 
#include "iceCompModel.hh"

PetscErrorCode IceCompModel::temperatureStep() {
  PetscErrorCode  ierr;

  if ((testname == 'F') || (testname == 'G')) {
    ierr = VecAXPY(vSigma, 1.0, vSigmaComp); CHKERRQ(ierr);   // Sigma = Sigma + Sigma_c
    ierr = IceModel::temperatureStep(); CHKERRQ(ierr);
    ierr = VecAXPY(vSigma, -1.0, vSigmaComp); CHKERRQ(ierr);  // Sigma = Sigma - Sigma_c
  } else {
    ierr = IceModel::temperatureStep(); CHKERRQ(ierr);
  }
  return 0;
}


// boundary conditions for tests F, G (same as EISMINT II Experiment F)
PetscScalar IceCompModel::Ggeo = 0.042;
PetscScalar IceCompModel::ST = 1.67e-5;
PetscScalar IceCompModel::Tmin = 223.15;  // K
PetscScalar IceCompModel::LforFG = 750000; // m
PetscScalar IceCompModel::ApforG = 200; // m


PetscErrorCode IceCompModel::createCompVecs() {
  PetscErrorCode ierr = DACreateLocalVector(grid.da3, &vSigmaComp); CHKERRQ(ierr);
  ierr = VecSet(vSigmaComp,0.0); CHKERRQ(ierr);
  compVecsCreated = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceCompModel::destroyCompVecs() {
  PetscErrorCode  ierr = VecDestroy(vSigmaComp); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::createCompViewers() {
  PetscErrorCode ierr;

  // must be called after IceModel::createViewers because diagnostic needs to be filled
  if ((testname=='F') || (testname=='G')) {
    ierr = createOneViewerIfDesired(&SigmaCompView, 'P',
                   "Sigma_C (comPensatory heat; K/a) at kd");  CHKERRQ(ierr);
  } else SigmaCompView = PETSC_NULL;
  
  // take over SigmaMapView to show only strain heating and not sum Sigma + Sigma_C
  if (runtimeViewers[cIndex('S')] != PETSC_NULL) {
    ierr = PetscViewerDestroy(runtimeViewers[cIndex('S')]); CHKERRQ(ierr);
    runtimeViewers[cIndex('S')] = PETSC_NULL;
    ierr = createOneViewerIfDesired(&compSigmaMapView, 'S',
                   "Sigma (strain heat; K/a) at kd");  CHKERRQ(ierr);
  } else compSigmaMapView = PETSC_NULL;
   
  compViewersCreated = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceCompModel::destroyCompViewers() {
  PetscErrorCode ierr;
  if (SigmaCompView != PETSC_NULL) { ierr = PetscViewerDestroy(SigmaCompView); CHKERRQ(ierr); }
  if (compSigmaMapView != PETSC_NULL) { ierr = PetscViewerDestroy(compSigmaMapView); CHKERRQ(ierr); }
  return 0;
}


PetscErrorCode IceCompModel::initTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **H, **accum, **Ts, ***T, ***Tb;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4;
  z=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  ierr = VecSet(vbed, 0); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  ierr = VecSet(vGhf, Ggeo); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      Ts[i][j] = Tmin + ST * r;
      if (r > LforFG - 1.0) { // if (essentially) outside of sheet
        H[i][j]=0.0;
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k=0; k<Mz; k++)
          T[i][j][k]=Ts[i][j];
      } else {
        r=PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
           bothexact(0.0,r,z,Mz,0.0,
                     &H[i][j],&accum[i][j],T[i][j],dummy1,dummy2,dummy3,dummy4);
        } else {
           bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                     &H[i][j],&accum[i][j],T[i][j],dummy1,dummy2,dummy3,dummy4);
        }
      }
      // fill with basal temp increased by geothermal flux
      for (PetscInt k=0; k<grid.p->Mbz; k++)
        Tb[i][j][k] = T[i][j][0] + bedrock.k * (grid.p->Mbz - k - 1) * grid.p->dz * Ggeo;
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);

  ierr = VecCopy(vH,vh); CHKERRQ(ierr);

  delete [] z;  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  return 0;
}


PetscErrorCode IceCompModel::getCompSourcesTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **accum, ***SigmaC;
  PetscScalar     dummy0;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4;
  z=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;

  // before temperature and flow step, set Sigma_c and accumulation from exact values
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k=0; k<Mz; k++)  SigmaC[i][j][k]=0.0;
      } else {
        r=PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
          bothexact(0.0,r,z,Mz,0.0,
                    &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
        } else {
          bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                    &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  delete [] z;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestFG() {
  // fills Vecs vH, vh, vAccum, vT, vu, vv, vw, vSigma, vSigmaComp
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **H, **accum, ***T, ***w, ***u, ***v, ***Sigma, ***SigmaC;
  PetscScalar     Ts, *z, *Uradial;

  z = new PetscScalar[Mz];
  Uradial = new PetscScalar[Mz];

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        H[i][j] = 0.0;
        Ts = Tmin + ST * r;
        for (PetscInt k=0; k<Mz; k++) {
          T[i][j][k] = Ts;
          u[i][j][k] = 0.0;
          v[i][j][k] = 0.0;
          w[i][j][k] = 0.0;
          Sigma[i][j][k] = 0.0;
          SigmaC[i][j][k] = 0.0;
        }
      } else {  // inside the sheet
        r = PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
          bothexact(0.0,r,z,Mz,0.0,
                    &H[i][j],&accum[i][j],T[i][j],Uradial,w[i][j],
                    Sigma[i][j],SigmaC[i][j]);
        } else {
          bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                    &H[i][j],&accum[i][j],T[i][j],Uradial,w[i][j],
                    Sigma[i][j],SigmaC[i][j]);
        }
        for (PetscInt k=0; k<Mz; k++) {
          u[i][j][k] = Uradial[k]*(xx/r);
          v[i][j][k] = Uradial[k]*(yy/r);
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  delete [] z;  delete [] Uradial;

  ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
  ierr = VecCopy(vH,vh); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceCompModel::updateCompViewers() {
  PetscErrorCode ierr;

  ierr = updateViewers();  CHKERRQ(ierr);

  Vec myg2;
  ierr = DACreateGlobalVector(grid.da2, &myg2); CHKERRQ(ierr);
  Vec* myvWork2d;
  ierr = VecDuplicateVecs(vh, 3, &myvWork2d); CHKERRQ(ierr);
  
  if (SigmaCompView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vSigmaComp, myg2, kd);
    ierr = VecScale(myg2, secpera); CHKERRQ(ierr);
    ierr = VecView(myg2, SigmaCompView); CHKERRQ(ierr);
  }  
  // now redraw Sigma view to be just the strain-heating, not the sum of strain-heating 
  //      and compensatory
  if (compSigmaMapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vSigma, myvWork2d[0], kd);
    ierr = getHorSliceOf3D(vSigmaComp, myvWork2d[1], kd);
    // Work2d[2] = Sigma - SigmaComp:
    ierr = VecWAXPY(myvWork2d[2],-1.0,myvWork2d[1],myvWork2d[0]); CHKERRQ(ierr);
    ierr = DALocalToGlobal(grid.da2, myvWork2d[2], INSERT_VALUES, myg2); CHKERRQ(ierr);
    ierr = VecScale(myg2, secpera); CHKERRQ(ierr);
    ierr = VecView(myg2, compSigmaMapView); CHKERRQ(ierr);
  }

  ierr = VecDestroyVecs(myvWork2d, 3); CHKERRQ(ierr);
  ierr = VecDestroy(myg2); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::computeTemperatureErrors(PetscScalar &gmaxTerr, PetscScalar &gavTerr) {

  PetscErrorCode ierr;
  PetscScalar    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  PetscScalar    **H, ***T;
  const PetscInt Mz = grid.p->Mz;
  
  PetscScalar   *z, *dummy1, *dummy2, *dummy3, *dummy4, *Tex;
  PetscScalar   junk0, junk1;
  z = new PetscScalar[Mz];  Tex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];
  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,z,Mz,0.0,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          default:  SETERRQ(1,"temperature errors only computable for tests F and G\n");
        }
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; k++) {  // only eval error if below num surface
          const PetscScalar Terr = PetscAbs(T[i][j][k] - Tex[k]);
          maxTerr = PetscMax(maxTerr,Terr);
          avcount += 1.0;
          avTerr += Terr;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);

  delete [] z;  delete [] Tex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PetscGlobalMax(&maxTerr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PetscGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeIceBedrockTemperatureErrors(
                                PetscScalar &gmaxTerr, PetscScalar &gavTerr,
                                PetscScalar &gmaxTberr, PetscScalar &gavTberr) {
  PetscErrorCode ierr;

  if (testname != 'K')
    SETERRQ(1,"ice and bedrock temperature errors only computable for test K\n");

  PetscScalar    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  PetscScalar    maxTberr = 0.0, avTberr = 0.0, avbcount = 0.0;
  PetscScalar    ***T, ***Tb, *Tex, *Tbex;
  const PetscInt    Mz = grid.p->Mz, Mbz = grid.p->Mbz;
  const PetscScalar dz = grid.p->dz;
 
  Tex = new PetscScalar[Mz];  Tbex = new PetscScalar[Mbz];

  for (PetscInt k=0; k<Mz; k++) {
    ierr = exactK(grid.p->year * secpera, k * dz, &Tex[k],(bedrock_is_ice_forK==PETSC_TRUE));
             CHKERRQ(ierr);
  }
  for (PetscInt k=0; k<Mbz; k++) {
    const PetscScalar depth = ((Mbz-1) - k) * dz;
    ierr = exactK(grid.p->year * secpera, -depth, &Tbex[k],(bedrock_is_ice_forK==PETSC_TRUE));
             CHKERRQ(ierr);
  }
    
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      for (PetscInt kb=0; kb<Mbz; kb++) { 
        const PetscScalar Tberr = PetscAbs(Tb[i][j][kb] - Tbex[kb]);
        maxTberr = PetscMax(maxTberr,Tberr);
        avbcount += 1.0;
        avTberr += Tberr;
      }
      for (PetscInt k=0; k<Mz; k++) { 
        const PetscScalar Terr = PetscAbs(T[i][j][k] - Tex[k]);
        maxTerr = PetscMax(maxTerr,Terr);
        avcount += 1.0;
        avTerr += Terr;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  delete [] Tex;  delete [] Tbex;
  
  ierr = PetscGlobalMax(&maxTerr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PetscGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero

  ierr = PetscGlobalMax(&maxTberr, &gmaxTberr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTberr, &gavTberr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavbcount;
  ierr = PetscGlobalSum(&avbcount, &gavbcount, grid.com); CHKERRQ(ierr);
  gavTberr = gavTberr/PetscMax(gavbcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeBasalTemperatureErrors(
      PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr) {

  PetscErrorCode  ierr;
  PetscScalar     ***T, domeT, domeTexact, Terr, avTerr;

  PetscScalar     dummy, z, Texact, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  domeT=0; domeTexact = 0; Terr=0; avTerr=0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      switch (testname) {
        case 'F':
          if (r > LforFG - 1.0) {  // outside of sheet
            Texact=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(0.0,r,&z,1,0.0,
                      &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        case 'G':
          if (r > LforFG -1.0) {  // outside of sheet
            Texact=Tmin + ST * r;  // = Ts
          } else {
            r=PetscMax(r,1.0);
            z=0.0;
            bothexact(grid.p->year*secpera,r,&z,1,ApforG,
                      &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        default:  SETERRQ(1,"temperature errors only computable for tests F and G\n");
      }

      if (i == (grid.p->Mx - 1)/2 && j == (grid.p->My - 1)/2) {
        domeT = T[i][j][0];
        domeTexact = Texact;
      }
      // compute maximum errors
      Terr = PetscMax(Terr,PetscAbsReal(T[i][j][0] - Texact));
      // add to sums for average errors
      avTerr += PetscAbs(T[i][j][0] - Texact);
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  
  PetscScalar gdomeT, gdomeTexact;

  ierr = PetscGlobalMax(&Terr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalMax(&domeT, &gdomeT, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&domeTexact, &gdomeTexact, grid.com); CHKERRQ(ierr);  
  centerTerr = PetscAbsReal(gdomeT - gdomeTexact);
  
  return 0;
}


PetscErrorCode IceCompModel::computeSigmaErrors(PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr) {

  PetscErrorCode ierr;
  PetscScalar    maxSigerr = 0.0, avSigerr = 0.0, avcount = 0.0;
  PetscScalar    **H, ***Sig;
  const PetscInt Mz = grid.p->Mz;
  
  PetscScalar   *z, *dummy1, *dummy2, *dummy3, *dummy4, *Sigex;
  PetscScalar   junk0, junk1;
  z = new PetscScalar[Mz];  Sigex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];
  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
    
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sig); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,z,Mz,0.0,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,dummy4);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,dummy4);
            break;
          default:  SETERRQ(1,"strain-heating (Sigma) errors only computable for tests F and G\n");
        }
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; k++) {  // only eval error if below num surface
          const PetscScalar Sigerr = PetscAbs(Sig[i][j][k] - Sigex[k]);
          maxSigerr = PetscMax(maxSigerr,Sigerr);
          avcount += 1.0;
          avSigerr += Sigerr;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sig); CHKERRQ(ierr);

  delete [] z;  delete [] Sigex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PetscGlobalMax(&maxSigerr, &gmaxSigmaerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avSigerr, &gavSigmaerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PetscGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavSigmaerr = gavSigmaerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeSurfaceVelocityErrors(
        PetscScalar &gmaxUerr, PetscScalar &gavUerr,
        PetscScalar &gmaxWerr, PetscScalar &gavWerr) {

  PetscErrorCode ierr;
  PetscScalar    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;
  PetscScalar    **H, ***unum, ***vnum, ***wnum;

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vu, &unum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vv, &vnum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vw, &wnum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        PetscScalar radialUex,wex;
        PetscScalar dummy0,dummy1,dummy2,dummy3,dummy4;
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        PetscScalar z = ks * grid.p->dz;
        switch (testname) {
          case 'F':
            bothexact(0.0,r,&z,1,0.0,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,&z,1,ApforG,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          default:  SETERRQ(1,"surface velocity errors only computed for tests F and G\n");
        }
        const PetscScalar uex = (xx/r) * radialUex;
        const PetscScalar vex = (yy/r) * radialUex;
        const PetscScalar Uerr = sqrt(PetscSqr(unum[i][j][ks] - uex)
                                      + PetscSqr(vnum[i][j][ks] - vex));
        maxUerr = PetscMax(maxUerr,Uerr);
        avUerr += Uerr;
        const PetscScalar Werr = PetscAbs(wnum[i][j][ks] - wex);
        maxWerr = PetscMax(maxWerr,Werr);
        avWerr += Werr;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vu, &unum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vv, &vnum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vw, &wnum); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.p->Mx*grid.p->My);
  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestK() {
  PetscErrorCode    ierr;
  const PetscInt    Mz = grid.p->Mz, Mbz = grid.p->Mbz; 
  const PetscScalar dz = grid.p->dz;
  PetscScalar       ***T, ***Tb;
  PetscScalar       *Tcol, *Tbcol;

  // evaluate exact solution in a column; all columns are the same
  Tcol = new PetscScalar[Mz];
  Tbcol = new PetscScalar[Mbz];
  for (PetscInt k=0; k<Mz; k++) {
    // evaluate and store in Tcol and Tbcol
    ierr = exactK(grid.p->year * secpera, k * dz, &Tcol[k],(bedrock_is_ice_forK==PETSC_TRUE));
             CHKERRQ(ierr);
  }
  for (PetscInt k=0; k<Mbz; k++) {
    // evaluate and store in Tcol and Tbcol
    const PetscScalar depth = ((Mbz-1) - k) * dz;
    ierr = exactK(grid.p->year * secpera, -depth, &Tbcol[k],(bedrock_is_ice_forK==PETSC_TRUE));
             CHKERRQ(ierr);
  }

  // copy column values into 3D arrays
  ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      for (PetscInt k=0; k<Mz; k++) {
        T[i][j][k] = Tcol[k];
      }
      for (PetscInt k=0; k<grid.p->Mbz; k++) {
        Tb[i][j][k] = Tbcol[k];
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3b, vTb, &Tb); CHKERRQ(ierr);

  delete [] Tcol;  delete [] Tbcol;

  // only communicate vT (as vTb will not be horizontally differentiated)
  ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::initTestK() {
  PetscErrorCode    ierr;

  ierr = VecSet(vbed, 0.0); CHKERRQ(ierr);
  ierr = VecSet(vAccum, 0.0); CHKERRQ(ierr);
  ierr = VecSet(vTs, 223.15); CHKERRQ(ierr);
  ierr = VecSet(vMask, MASK_SHEET); CHKERRQ(ierr);
  ierr = VecSet(vGhf, 0.042); CHKERRQ(ierr);
  ierr = VecSet(vH, 3000.0); CHKERRQ(ierr);
  ierr = VecSet(vHmelt, 0.0); CHKERRQ(ierr);
  ierr = VecSet(vtau, 0.0); CHKERRQ(ierr);
  ierr = VecCopy(vH,vh); CHKERRQ(ierr);

  ierr = fillSolnTestK(); CHKERRQ(ierr);
  return 0;
}

