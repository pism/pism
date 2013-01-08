// Copyright (C) 2004-2013 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <petscdmda.h>
#include "tests/exactTestsFG.h" 
#include "tests/exactTestK.h" 
#include "tests/exactTestO.h" 
#include "iceCompModel.hh"
#include "PISMStressBalance.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

// boundary conditions for tests F, G (same as EISMINT II Experiment F)
const PetscScalar IceCompModel::Ggeo = 0.042;
const PetscScalar IceCompModel::ST = 1.67e-5;
const PetscScalar IceCompModel::Tmin = 223.15;  // K
const PetscScalar IceCompModel::LforFG = 750000; // m
const PetscScalar IceCompModel::ApforG = 200; // m


/*! Re-implemented so that we can add compensatory Sigma in Tests F and G. */
PetscErrorCode IceCompModel::temperatureStep(
       PetscScalar* vertSacrCount, PetscScalar* bulgeCount) {
  PetscErrorCode  ierr;

  if ((testname == 'F') || (testname == 'G')) {
    IceModelVec3 *Sigma3;
    ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr);

    ierr = Sigma3->add(1.0, SigmaComp3); CHKERRQ(ierr);	// Sigma = Sigma + Sigma_c
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount); CHKERRQ(ierr);
    ierr = Sigma3->add(-1.0, SigmaComp3); CHKERRQ(ierr); // Sigma = Sigma - Sigma_c
  } else {
    ierr = IceModel::temperatureStep(vertSacrCount,bulgeCount); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceCompModel::initTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.Mz;
  PetscScalar     **H, **accum, **Ts;
  PetscScalar     *dummy1, *dummy2, *dummy3, *dummy4;

  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  ierr = vbed.set(0); CHKERRQ(ierr);
  ierr = vGhf.set(Ggeo); CHKERRQ(ierr);

  PetscScalar *T;
  T = new PetscScalar[grid.Mz];

  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = artm.get_array(Ts); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r = grid.radius(i,j);
      Ts[i][j] = Tmin + ST * r;
      if (r > LforFG - 1.0) { // if (essentially) outside of sheet
        H[i][j] = 0.0;
        accum[i][j] = -ablationRateOutside/secpera;
        for (PetscInt k = 0; k < Mz; k++)
          T[k]=Ts[i][j];
      } else {
        r = PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
          bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                     &H[i][j],&accum[i][j],T,dummy1,dummy2,dummy3,dummy4);
        } else {
          bothexact(grid.time->current(),r,&grid.zlevels[0],Mz,ApforG,
                     &H[i][j],&accum[i][j],T,dummy1,dummy2,dummy3,dummy4);
        }
      }
      ierr = T3.setInternalColumn(i,j,T); CHKERRQ(ierr);

    }
  }

  ierr =     vH.end_access(); CHKERRQ(ierr);
  ierr =     T3.end_access(); CHKERRQ(ierr);

  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);

  ierr = vH.update_ghosts(); CHKERRQ(ierr);
  
  ierr = T3.update_ghosts(); CHKERRQ(ierr);

  ierr = vH.copy_to(vh); CHKERRQ(ierr);

  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  delete [] T;
  
  return 0;
}


PetscErrorCode IceCompModel::getCompSourcesTestFG() {
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.Mz;
  PetscScalar     **accum;
  PetscScalar     dummy0;
  PetscScalar     *dummy1, *dummy2, *dummy3, *dummy4;

  dummy1=new PetscScalar[Mz];  dummy2=new PetscScalar[Mz];
  dummy3=new PetscScalar[Mz];  dummy4=new PetscScalar[Mz];

  PetscScalar *SigmaC;
  SigmaC = new PetscScalar[Mz];

  const PetscScalar
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  // before temperature and flow step, set Sigma_c and accumulation from exact values
  ierr = acab.get_array(accum); CHKERRQ(ierr);
  ierr = SigmaComp3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r = grid.radius(i,j);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        ierr = SigmaComp3.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else {
        r = PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
          bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                    &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC);
        } else {
          bothexact(grid.time->current(),r,&grid.zlevels[0],Mz,ApforG,
                    &dummy0,&accum[i][j],dummy1,dummy2,dummy3,dummy4,SigmaC);
        }
        for (PetscInt k=0;  k<Mz;  k++) // scale Sigma to J/(s m^3)
          SigmaC[k] = SigmaC[k] * ice_rho * ice_c;
        ierr = SigmaComp3.setInternalColumn(i,j,SigmaC); CHKERRQ(ierr);
      }
    }
  }
  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = SigmaComp3.end_access(); CHKERRQ(ierr);

  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  delete [] SigmaC;

  return 0;
}


PetscErrorCode IceCompModel::fillSolnTestFG() {
  // fills Vecs vH, vh, vAccum, T3, u3, v3, w3, Sigma3, vSigmaComp
  PetscErrorCode  ierr;
  PetscInt        Mz=grid.Mz;
  PetscScalar     **H, **accum;
  PetscScalar     Ts, *Uradial;

  IceModelVec3 *u3, *v3, *w3, *Sigma3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);
  ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr);

  Uradial = new PetscScalar[Mz];

  PetscScalar *T, *u, *v, *w, *Sigma, *SigmaC;
  T = new PetscScalar[grid.Mz];
  u = new PetscScalar[grid.Mz];
  v = new PetscScalar[grid.Mz];
  w = new PetscScalar[grid.Mz];
  Sigma = new PetscScalar[grid.Mz];
  SigmaC = new PetscScalar[grid.Mz];

  const PetscScalar
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = acab.get_array(accum); CHKERRQ(ierr);
  
  ierr = T3.begin_access(); CHKERRQ(ierr);

  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = Sigma3->begin_access(); CHKERRQ(ierr);
  ierr = SigmaComp3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar xx = grid.x[i], yy = grid.y[j], r = grid.radius(i,j);
      if (r > LforFG - 1.0) {  // outside of sheet
        accum[i][j] = -ablationRateOutside/secpera;
        H[i][j] = 0.0;
        Ts = Tmin + ST * r;
        ierr = T3.setColumn(i,j,Ts); CHKERRQ(ierr);
        ierr = u3->setColumn(i,j,0.0); CHKERRQ(ierr);
        ierr = v3->setColumn(i,j,0.0); CHKERRQ(ierr);
        ierr = w3->setColumn(i,j,0.0); CHKERRQ(ierr);
        ierr = Sigma3->setColumn(i,j,0.0); CHKERRQ(ierr);
        ierr = SigmaComp3.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else {  // inside the sheet
        r = PetscMax(r,1.0); // avoid singularity at origin
        if (testname == 'F') {
          bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                    &H[i][j],&accum[i][j],T,Uradial,w, Sigma,SigmaC);
        } else {
          bothexact(grid.time->current(),r,&grid.zlevels[0],Mz,ApforG,
                    &H[i][j],&accum[i][j],T,Uradial,w, Sigma,SigmaC);
        }
        for (PetscInt k = 0; k < Mz; k++) {
          u[k] = Uradial[k]*(xx/r);
          v[k] = Uradial[k]*(yy/r);
          Sigma[k] = Sigma[k] * ice_rho * ice_c; // scale Sigma to J/(s m^3)
          SigmaC[k] = SigmaC[k] * ice_rho * ice_c; // scale SigmaC to J/(s m^3)
        }
        ierr = T3.setInternalColumn(i,j,T); CHKERRQ(ierr);
        ierr = u3->setInternalColumn(i,j,u); CHKERRQ(ierr);
        ierr = v3->setInternalColumn(i,j,v); CHKERRQ(ierr);
        ierr = w3->setInternalColumn(i,j,w); CHKERRQ(ierr);
        ierr = Sigma3->setInternalColumn(i,j,Sigma); CHKERRQ(ierr);
        ierr = SigmaComp3.setInternalColumn(i,j,SigmaC); CHKERRQ(ierr);
      }      
    }
  }

  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);
  ierr = SigmaComp3.end_access(); CHKERRQ(ierr);
  
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = acab.end_access(); CHKERRQ(ierr);

  delete [] Uradial;

  delete [] T;  delete [] u;  delete [] v;  delete [] w;
  delete [] Sigma;  delete [] SigmaC;

  ierr = vH.update_ghosts(); CHKERRQ(ierr);
  ierr = vH.copy_to(vh); CHKERRQ(ierr);

  ierr = T3.update_ghosts(); CHKERRQ(ierr);

  ierr = u3->update_ghosts(); CHKERRQ(ierr);

  ierr = v3->update_ghosts(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceCompModel::computeTemperatureErrors(
                                  PetscScalar &gmaxTerr, PetscScalar &gavTerr) {

  PetscErrorCode ierr;
  PetscScalar    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  PetscScalar    **H;
  const PetscInt Mz = grid.Mz;
  
  PetscScalar   *dummy1, *dummy2, *dummy3, *dummy4, *Tex;
  PetscScalar   junk0, junk1;
  
  Tex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];

  PetscScalar *T;

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r = grid.radius(i,j);
      ierr = T3.getInternalColumn(i,j,&T); CHKERRQ(ierr);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          case 'G':
            bothexact(grid.time->current(),r,&grid.zlevels[0],Mz,ApforG,
                      &junk0,&junk1,Tex,dummy1,dummy2,dummy3,dummy4);
            break;
          default:  SETERRQ(grid.com, 1,"temperature errors only computable for tests F and G\n");
        }
        const PetscInt ks = grid.kBelowHeight(H[i][j]);
        for (PetscInt k = 0; k < ks; k++) {  // only eval error if below num surface
          const PetscScalar Terr = PetscAbs(T[k] - Tex[k]);
          maxTerr = PetscMax(maxTerr,Terr);
          avcount += 1.0;
          avTerr += Terr;
        }
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = T3.end_access(); CHKERRQ(ierr);

  delete [] Tex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PISMGlobalMax(&maxTerr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PISMGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeIceBedrockTemperatureErrors(
                                PetscScalar &gmaxTerr, PetscScalar &gavTerr,
                                PetscScalar &gmaxTberr, PetscScalar &gavTberr) {
  PetscErrorCode ierr;

  if ((testname != 'K') && (testname != 'O'))
    SETERRQ(grid.com, 1,"ice and bedrock temperature errors only computable for tests K and O\n");

  PetscScalar    maxTerr = 0.0, avTerr = 0.0, avcount = 0.0;
  PetscScalar    maxTberr = 0.0, avTberr = 0.0, avbcount = 0.0;
  const PetscInt    Mz = grid.Mz;
 
  PetscScalar    *Tex, *Tbex, *T, *Tb;
  PetscScalar    FF;
  Tex = new PetscScalar[Mz];  

  IceModelVec3BTU *bedrock_temp;

  BTU_Verification *my_btu = dynamic_cast<BTU_Verification*>(btu);
  if (my_btu == NULL) SETERRQ(grid.com, 1, "my_btu == NULL");
  ierr = my_btu->get_temp(bedrock_temp); CHKERRQ(ierr);

  vector<double> zblevels = bedrock_temp->get_levels();
  int Mbz = (int)zblevels.size();
  Tbex = new PetscScalar[Mbz];

  switch (testname) {
    case 'K':
      for (PetscInt k = 0; k < Mz; k++) {
        ierr = exactK(grid.time->current(), grid.zlevels[k], &Tex[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      for (PetscInt k = 0; k < Mbz; k++) {
        ierr = exactK(grid.time->current(), zblevels[k], &Tbex[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      break;
    case 'O':
      PetscScalar dum1, dum2, dum3, dum4;
      for (PetscInt k = 0; k < Mz; k++) {
        ierr = exactO(grid.zlevels[k], &Tex[k], &dum1, &dum2, &dum3, &dum4);
             CHKERRQ(ierr);
      }
      for (PetscInt k = 0; k < Mbz; k++) {
        ierr = exactO(zblevels[k], &Tbex[k], &dum1, &dum2, &dum3, &dum4);
             CHKERRQ(ierr);
      }
      break;
    default: SETERRQ(grid.com, 2,"again: ice and bedrock temperature errors only for tests K and O\n");
  }
    
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = bedrock_temp->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      ierr = bedrock_temp->getInternalColumn(i,j,&Tb); CHKERRQ(ierr);
      for (PetscInt kb = 0; kb < Mbz; kb++) { 
        const PetscScalar Tberr = PetscAbs(Tb[kb] - Tbex[kb]);
        maxTberr = PetscMax(maxTberr,Tberr);
        avbcount += 1.0;
        avTberr += Tberr;
      }
      ierr = T3.getInternalColumn(i,j,&T); CHKERRQ(ierr);
      for (PetscInt k = 0; k < Mz; k++) { 
        const PetscScalar Terr = PetscAbs(T[k] - Tex[k]);
        maxTerr = PetscMax(maxTerr,Terr);
        avcount += 1.0;
        avTerr += Terr;
      }
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = bedrock_temp->end_access(); CHKERRQ(ierr);

  delete [] Tex;  delete [] Tbex;
  
  ierr = PISMGlobalMax(&maxTerr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PISMGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/PetscMax(gavcount,1.0);  // avoid div by zero

  ierr = PISMGlobalMax(&maxTberr, &gmaxTberr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avTberr, &gavTberr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavbcount;
  ierr = PISMGlobalSum(&avbcount, &gavbcount, grid.com); CHKERRQ(ierr);
  gavTberr = gavTberr/PetscMax(gavbcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeBasalTemperatureErrors(
      PetscScalar &gmaxTerr, PetscScalar &gavTerr, PetscScalar &centerTerr) {

  PetscErrorCode  ierr;
  PetscScalar     domeT, domeTexact, Terr, avTerr;

  PetscScalar     dummy, z, Texact, dummy1, dummy2, dummy3, dummy4, dummy5;

  ierr = T3.begin_access(); CHKERRQ(ierr);

  domeT=0; domeTexact = 0; Terr=0; avTerr=0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar r = grid.radius(i,j);
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
            bothexact(grid.time->current(),r,&z,1,ApforG,
                      &dummy5,&dummy,&Texact,&dummy1,&dummy2,&dummy3,&dummy4);
          }
          break;
        default:  SETERRQ(grid.com, 1,"temperature errors only computable for tests F and G\n");
      }

      const PetscScalar Tbase = T3.getValZ(i,j,0.0);
      if (i == (grid.Mx - 1)/2 && j == (grid.My - 1)/2) {
        domeT = Tbase;
        domeTexact = Texact;
      }
      // compute maximum errors
      Terr = PetscMax(Terr,PetscAbsReal(Tbase - Texact));
      // add to sums for average errors
      avTerr += PetscAbs(Tbase - Texact);
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  
  PetscScalar gdomeT, gdomeTexact;

  ierr = PISMGlobalMax(&Terr, &gmaxTerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avTerr, &gavTerr, grid.com); CHKERRQ(ierr);
  gavTerr = gavTerr/(grid.Mx*grid.My);
  ierr = PISMGlobalMax(&domeT, &gdomeT, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&domeTexact, &gdomeTexact, grid.com); CHKERRQ(ierr);  
  centerTerr = PetscAbsReal(gdomeT - gdomeTexact);
  
  return 0;
}


PetscErrorCode IceCompModel::computeSigmaErrors(
                      PetscScalar &gmaxSigmaerr, PetscScalar &gavSigmaerr) {

  PetscErrorCode ierr;
  PetscScalar    maxSigerr = 0.0, avSigerr = 0.0, avcount = 0.0;
  PetscScalar    **H;
  const PetscInt Mz = grid.Mz;
  
  PetscScalar   *dummy1, *dummy2, *dummy3, *dummy4, *Sigex;
  PetscScalar   junk0, junk1;
  
  Sigex = new PetscScalar[Mz];
  dummy1 = new PetscScalar[Mz];  dummy2 = new PetscScalar[Mz];
  dummy3 = new PetscScalar[Mz];  dummy4 = new PetscScalar[Mz];

  const PetscScalar
    ice_rho   = config.get("ice_density"),
    ice_c     = config.get("ice_specific_heat_capacity");

  PetscScalar *Sig;
  IceModelVec3 *Sigma3;
  ierr = stress_balance->get_volumetric_strain_heating(Sigma3); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = Sigma3->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar r = grid.radius(i,j);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        switch (testname) {
          case 'F':
            bothexact(0.0,r,&grid.zlevels[0],Mz,0.0,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,dummy4);
            break;
          case 'G':
            bothexact(grid.time->current(),r,&grid.zlevels[0],Mz,ApforG,
                      &junk0,&junk1,dummy1,dummy2,dummy3,Sigex,dummy4);
            break;
          default:
            SETERRQ(grid.com, 1,"strain-heating (Sigma) errors only computable for tests F and G\n");
        }
        for (PetscInt k = 0; k < Mz; k++)
          Sigex[k] = Sigex[k] * ice_rho * ice_c; // scale exact Sigma to J/(s m^3)
        const PetscInt ks = grid.kBelowHeight(H[i][j]);
        ierr = Sigma3->getInternalColumn(i,j,&Sig); CHKERRQ(ierr);
        for (PetscInt k = 0; k < ks; k++) {  // only eval error if below num surface
          const PetscScalar Sigerr = PetscAbs(Sig[k] - Sigex[k]);
          maxSigerr = PetscMax(maxSigerr,Sigerr);
          avcount += 1.0;
          avSigerr += Sigerr;
        }
      }
    }
  }
  ierr =     vH.end_access(); CHKERRQ(ierr);
  ierr = Sigma3->end_access(); CHKERRQ(ierr);

  delete [] Sigex;
  delete [] dummy1;  delete [] dummy2;  delete [] dummy3;  delete [] dummy4;
  
  ierr = PISMGlobalMax(&maxSigerr, &gmaxSigmaerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avSigerr, &gavSigmaerr, grid.com); CHKERRQ(ierr);
  PetscScalar  gavcount;
  ierr = PISMGlobalSum(&avcount, &gavcount, grid.com); CHKERRQ(ierr);
  gavSigmaerr = gavSigmaerr/PetscMax(gavcount,1.0);  // avoid div by zero
  return 0;
}


PetscErrorCode IceCompModel::computeSurfaceVelocityErrors(
        PetscScalar &gmaxUerr, PetscScalar &gavUerr,
        PetscScalar &gmaxWerr, PetscScalar &gavWerr) {

  PetscErrorCode ierr;
  PetscScalar    maxUerr = 0.0, maxWerr = 0.0, avUerr = 0.0, avWerr = 0.0;
  PetscScalar    **H;

  IceModelVec3 *u3, *v3, *w3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      PetscScalar xx = grid.x[i], yy = grid.y[j], r = grid.radius(i,j);
      if ((r >= 1.0) && (r <= LforFG - 1.0)) {  // only evaluate error if inside sheet 
                                                // and not at central singularity
        PetscScalar radialUex,wex;
        PetscScalar dummy0,dummy1,dummy2,dummy3,dummy4;
        switch (testname) {
          case 'F':
            bothexact(0.0,r,&(H[i][j]),1,0.0,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          case 'G':
            bothexact(grid.time->current(),r,&(H[i][j]),1,ApforG,
                      &dummy0,&dummy1,&dummy2,&radialUex,&wex,&dummy3,&dummy4);
            break;
          default:  SETERRQ(grid.com, 1,"surface velocity errors only computed for tests F and G\n");
        }
        const PetscScalar uex = (xx/r) * radialUex;
        const PetscScalar vex = (yy/r) * radialUex;
        // note that because getValZ does linear interpolation and H[i][j] is not exactly at
        // a grid point, this causes nonzero errors even with option -eo
        const PetscScalar Uerr = sqrt(PetscSqr(u3->getValZ(i,j,H[i][j]) - uex)
                                      + PetscSqr(v3->getValZ(i,j,H[i][j]) - vex));
        maxUerr = PetscMax(maxUerr,Uerr);
        avUerr += Uerr;
        const PetscScalar Werr = PetscAbs(w3->getValZ(i,j,H[i][j]) - wex);
        maxWerr = PetscMax(maxWerr,Werr);
        avWerr += Werr;
      }
    }
  }
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxUerr, &gmaxUerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&maxWerr, &gmaxWerr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&avUerr, &gavUerr, grid.com); CHKERRQ(ierr);
  gavUerr = gavUerr/(grid.Mx*grid.My);
  ierr = PISMGlobalSum(&avWerr, &gavWerr, grid.com); CHKERRQ(ierr);
  gavWerr = gavWerr/(grid.Mx*grid.My);
  return 0;
}


PetscErrorCode IceCompModel::computeBasalMeltRateErrors(
                   PetscScalar &gmaxbmelterr, PetscScalar &gminbmelterr) {
  PetscErrorCode ierr;
  PetscScalar    maxbmelterr = -9.99e40, minbmelterr = 9.99e40, err;
  PetscScalar    bmelt,dum1,dum2,dum3,dum4;

  if (testname != 'O')
    SETERRQ(grid.com, 1,"basal melt rate errors are only computable for test O\n");

  // we just need one constant from exact solution:
  ierr = exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt); CHKERRQ(ierr);

  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; i++) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; j++) {
      err = PetscAbs(vbmr(i,j) - bmelt);
      maxbmelterr = PetscMax(maxbmelterr, err);
      minbmelterr = PetscMin(minbmelterr, err);
    }
  }
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&maxbmelterr, &gmaxbmelterr, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMin(&minbmelterr, &gminbmelterr, grid.com); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillTemperatureSolnTestsKO() {
  PetscErrorCode    ierr;
  const PetscInt    Mz = grid.Mz;

  PetscScalar       *Tcol;
  PetscScalar       dum1, dum2, dum3, dum4;
  PetscScalar    FF;
  Tcol = new PetscScalar[Mz];

  // evaluate exact solution in a column; all columns are the same
  switch (testname) {
    case 'K':
      for (PetscInt k=0; k<Mz; k++) {
        ierr = exactK(grid.time->current(), grid.zlevels[k], &Tcol[k], &FF,
                      (bedrock_is_ice_forK==PETSC_TRUE)); CHKERRQ(ierr);
      }
      break;
    case 'O':
      for (PetscInt k=0; k<Mz; k++) {
        ierr = exactO(grid.zlevels[k], &Tcol[k], &dum1, &dum2, &dum3, &dum4); CHKERRQ(ierr);
      }
      break;
    default:  SETERRQ(grid.com, 2,"only fills temperature solutions for tests K and O\n");
  }

  // copy column values into 3D arrays
  ierr = T3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.setInternalColumn(i,j,Tcol); CHKERRQ(ierr);
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);

  delete [] Tcol;

  // communicate T
  ierr = T3.update_ghosts(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::fillBasalMeltRateSolnTestO() {
  PetscErrorCode    ierr;
  PetscScalar       bmelt, dum1, dum2, dum3, dum4;
  if (testname != 'O') { SETERRQ(grid.com, 1,"only fills basal melt rate soln for test O\n"); }

  // we just need one constant from exact solution:
  ierr = exactO(0.0, &dum1, &dum2, &dum3, &dum4, &bmelt); CHKERRQ(ierr);

  ierr = vbmr.set(bmelt); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceCompModel::initTestsKO() {
  PetscErrorCode    ierr;

  if (testname == 'K') {
    bool Mbz_set;
    int Mbz;
    ierr = PISMOptionsInt("-Mbz", "Number of levels in the bedrock thermal model",
                          Mbz, Mbz_set); CHKERRQ(ierr);
    if (Mbz_set && Mbz < 2) {
      PetscPrintf(grid.com, "PISM ERROR: pismv test K requires a bedrock thermal layer 1000m deep.\n");
      PISMEnd();
    }
  }

  ierr = acab.set(0.0); CHKERRQ(ierr);
  ierr = artm.set(223.15); CHKERRQ(ierr);

  ierr = vbed.set(0.0); CHKERRQ(ierr);
  ierr = vGhf.set(0.042); CHKERRQ(ierr);
  ierr = vH.set(3000.0); CHKERRQ(ierr);
  ierr = vH.copy_to(vh); CHKERRQ(ierr);

  ierr = fillTemperatureSolnTestsKO(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BTU_Verification::get_temp(IceModelVec3BTU* &result) {
  result = &temp;
  return 0;
}

PetscErrorCode BTU_Verification::bootstrap() {
  PetscErrorCode ierr;

  if (Mbz < 2) return 0;

  vector<double> Tbcol(Mbz),
    zlevels = temp.get_levels();
  PetscScalar       dum1, dum2, dum3, dum4;
  PetscScalar    FF;

  // evaluate exact solution in a column; all columns are the same
  switch (testname) {
    case 'K':
      for (PetscInt k=0; k<Mbz; k++) {
        if (exactK(grid.time->current(), zlevels[k], &Tbcol[k], &FF,
                   (bedrock_is_ice==PETSC_TRUE)))
          SETERRQ1(grid.com, 1,"exactK() reports that level %9.7f is below B0 = -1000.0 m\n",
                   zlevels[k]);
      }
      break;
    case 'O':
      for (PetscInt k=0; k<Mbz; k++) {
        ierr = exactO(zlevels[k], &Tbcol[k], &dum1, &dum2, &dum3, &dum4); CHKERRQ(ierr);
      }
      break;
    default:
      {
        ierr = PISMBedThermalUnit::bootstrap(); CHKERRQ(ierr);
      }
  }

  // copy column values into 3D arrays
  ierr = temp.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = temp.setInternalColumn(i,j,&Tbcol[0]); CHKERRQ(ierr);
    }
  }
  ierr = temp.end_access(); CHKERRQ(ierr);

  return 0;
}
