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
#include "exact/exactTestsFG.h" 
#include "iceCompModel.hh"

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
PetscScalar IceCompModel::Ggeo = 0.042;
PetscScalar IceCompModel::ST = 1.67e-5;
PetscScalar IceCompModel::Tmin = 223.15;  // K
PetscScalar IceCompModel::LforFG = 750000; // m
PetscScalar IceCompModel::ApforG = 200; // m


PetscErrorCode IceCompModel::createCompVecs() {
  PetscErrorCode ierr = DACreateLocalVector(grid.da3, &vSigmaComp); CHKERRQ(ierr);
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
    ierr = createOneViewerIfDesired(&SigmaCompView, 'P',"Sigma_C (comPensatory heat; K/a) at kd");  CHKERRQ(ierr);
  } else SigmaCompView = PETSC_NULL;
  
  // take over SigmaMapView to show only strain heating and not sum Sigma + Sigma_C
  if (SigmaMapView != PETSC_NULL) {
    ierr = PetscViewerDestroy(SigmaMapView); CHKERRQ(ierr);
    SigmaMapView = PETSC_NULL;
    ierr = createOneViewerIfDesired(&compSigmaMapView, 'S',"Sigma (strain heat; K/a) at kd");  CHKERRQ(ierr);
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


PetscErrorCode IceCompModel::updateTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.p->Mz;
  PetscScalar     **H, **accum, ***T, ***w, ***u, ***v, ***Sigma, ***SigmaC;
  PetscScalar     Ts, dummy0;
  PetscScalar     *z, *dummy1, *dummy2, *dummy3, *dummy4, *Uradial;
  z=new PetscScalar[Mz];  Uradial=new PetscScalar[Mz];
  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  // before temperature and flow step, set Sigma_c and accumulation from
  // exact values for test F; also add SigmaC to Sigma for temperature computation
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vu, &u); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vv, &v); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da3, vw, &w); CHKERRQ(ierr);
  }

  for (PetscInt k=0; k<Mz; k++)
    z[k]=k*grid.p->dz;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k=0; k<Mz; k++)  SigmaC[i][j][k]=0.0;
        if (exactOnly == PETSC_TRUE) {
          H[i][j] = 0.0;
          Ts = Tmin + ST * r;
          for (PetscInt k=0; k<Mz; k++) {
            T[i][j][k]=Ts;
            u[i][j][k]=0.0;
            v[i][j][k]=0.0;
            w[i][j][k]=0.0;
            Sigma[i][j][k]=0.0;
          }
        }
      } else {
        r=PetscMax(r,1.0); // avoid singularity at origin
        if (exactOnly==PETSC_TRUE) {
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
            u[i][j][k]=Uradial[k]*(xx/r);
            v[i][j][k]=Uradial[k]*(yy/r);
          }
        } else { // actually do update with compensatory Sigma; not exactOnly
          if (testname == 'F') {
            bothexact(0.0,r,z,Mz,0.0,
                      &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
          } else {
             bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC[i][j]);
          }
          for (PetscInt k=0; k<Mz; k++)
            Sigma[i][j][k] += SigmaC[i][j][k]; // for rest of model, Sigma is both terms
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sigma); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigmaComp, &SigmaC); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigma, INSERT_VALUES, vSigma); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da3, vSigmaComp, INSERT_VALUES, vSigmaComp); CHKERRQ(ierr);

  if (exactOnly==PETSC_TRUE) {
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vT, &T); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vu, &u); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vv, &v); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da3, vw, &w); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da2, vH, INSERT_VALUES, vH); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vT, INSERT_VALUES, vT); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vu, INSERT_VALUES, vu); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vv, INSERT_VALUES, vv); CHKERRQ(ierr);
    ierr = DALocalToLocalBegin(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
    ierr = DALocalToLocalEnd(grid.da3, vw, INSERT_VALUES, vw); CHKERRQ(ierr);
    ierr = VecCopy(vH,vh); CHKERRQ(ierr);
  }

  delete [] z;  delete [] Uradial;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  return 0;
}


PetscErrorCode IceCompModel::updateCompViewers() {
  PetscErrorCode ierr;

  ierr = updateViewers();  CHKERRQ(ierr);
  if (SigmaCompView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vSigmaComp, g2, kd);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, SigmaCompView); CHKERRQ(ierr);
  }  
  // now redraw Sigma view to be just the strain-heating, not the sum of strain-heating 
  //      and compensatory
  if (compSigmaMapView != PETSC_NULL) {
    ierr = getHorSliceOf3D(vSigma, vWork2d[0], kd);
    ierr = getHorSliceOf3D(vSigmaComp, vWork2d[1], kd);
    ierr = VecWAXPY(vWork2d[2],-1.0,vWork2d[1],vWork2d[0]); CHKERRQ(ierr);  // Work2d[2] = Sigma - SigmaComp
    ierr = DALocalToGlobal(grid.da2, vWork2d[2], INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr = VecScale(g2, secpera); CHKERRQ(ierr);
    ierr = VecView(g2, compSigmaMapView); CHKERRQ(ierr);
  }
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
  
  PetscScalar   *z, *dummy1, *dummy2, *dummy3, *Sigex, *SigCompex;
  PetscScalar   junk0, junk1;
  z = new PetscScalar[Mz];  Sigex = new PetscScalar[Mz];  SigCompex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];
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
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,SigCompex);
            break;
          case 'G':
            bothexact(grid.p->year*secpera,r,z,Mz,ApforG,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,SigCompex);
            break;
          default:  SETERRQ(1,"strain-heating (Sigma) errors only computable for tests F and G\n");
        }
        // verbPrintf(1,grid.com,"%e|",Sigex[3]);
        const PetscInt ks = static_cast<PetscInt>(floor(H[i][j]/grid.p->dz));
        for (PetscInt k=0; k<ks; k++) {  // only eval error if below num surface
          const PetscScalar actualSignum  = Sig[i][j][k] - SigCompex[k];
          const PetscScalar Sigerr = PetscAbs(actualSignum - Sigex[k]);
          maxSigerr = PetscMax(maxSigerr,Sigerr);
          avcount += 1.0;
          avSigerr += Sigerr;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da3, vSigma, &Sig); CHKERRQ(ierr);

  delete [] z;  delete [] Sigex;  delete [] SigCompex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;
  
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
  PetscScalar    **H, **unum, **vnum, **wnum;

  ierr = getSurfaceValuesOf3D(vu, vWork2d[0]); CHKERRQ(ierr); // = numerical surface val of u
  ierr = getSurfaceValuesOf3D(vv, vWork2d[1]); CHKERRQ(ierr); // = numerical surface val of v
  ierr = getSurfaceValuesOf3D(vw, vWork2d[2]); CHKERRQ(ierr); // = numerical surface val of w

  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &unum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[1], &vnum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[2], &wnum); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r,xx,yy;
      mapcoords(i,j,xx,yy,r);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                               // and not at central singularity
        PetscScalar radialUex,wex;
        PetscScalar dummy0,dummy1,dummy2,dummy3,dummy4;
        PetscScalar z = H[i][j];
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
        // verbPrintf(1,grid.com,"[%f|%f]",radialUex,wex);
        const PetscScalar uex = (xx/r) * radialUex;
        const PetscScalar vex = (yy/r) * radialUex;
        const PetscScalar Uerr = sqrt(PetscSqr(unum[i][j] - uex) + PetscSqr(vnum[i][j] - vex));
        maxUerr = PetscMax(maxUerr,Uerr);
        avUerr += Uerr;
        const PetscScalar Werr = PetscAbs(wnum[i][j] - wex);
        maxWerr = PetscMax(maxWerr,Werr);
        avWerr += Werr;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &unum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[1], &vnum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[2], &wnum); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.p->Mx*grid.p->My);
  ierr = PetscGlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.p->Mx*grid.p->My);
  return 0;
}

