// Copyright (C) 2007-2008 Ed Bueler
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

#include <cstring>
#include <petsc.h>
#include "icePSTexModel.hh"

const int Nexpers = 6;

const int NAME_LENGTH = 10;
char exper_names[Nexpers][NAME_LENGTH] = { "P0A",
                                           "P0I",
                                           "P1",
                                           "P2",
                                           "P3",
                                           "P4" };

PetscScalar stream_angle[3] = {0.0, 100.0, 225.0};  // degrees

const PetscScalar DEFAULT_PHI_STRONG = 20.0;

PetscScalar exper_params[Nexpers][4][3] = {

  {// P0A
    {0.0, 0.0, 0.0},    // bed depth at end
    {70.0, 70.0, 70.0}, // stream_width in km; INACTIVE
    {5.0, 5.0, 5.0},    // upstream till phi; INACTIVE
    {5.0, 5.0, 5.0}     // downstream till phi; INACTIVE
  },

  {// P0I
    {1000.0, 1000.0, 1000.0},    // bed depth at end
    {70.0, 70.0, 70.0}, // stream_width in km; INACTIVE
    {5.0, 5.0, 5.0},    // upstream till phi; INACTIVE
    {5.0, 5.0, 5.0}     // downstream till phi; INACTIVE
  },

  {// P1
    {0.0, 0.0, 0.0},    // bed depth at end
    {70.0, 70.0, 70.0}, // stream_width in km
    {5.0, 5.0, 5.0},    // upstream till phi
    {5.0, 5.0, 5.0}     // downstream till phi
  },

  {// P2
    {0.0, 0.0, 0.0},     // bed depth at end
    {70.0, 40.0, 100.0}, // stream_width in km
    {5.0, 5.0, 5.0},     // upstream till phi
    {5.0, 5.0, 5.0}      // downstream till phi
  },

  {// P3: default
    {1000.0, 1000.0, 1000.0},  // bed depth at end
    {70.0, 40.0, 100.0},       // stream_width in km
    {5.0, 5.0, 5.0},           // upstream till phi
    {5.0, 5.0, 5.0}            // downstream till phi
  },

  {// P4
    {1000.0, 1000.0, 1000.0},  // bed depth at end
    {70.0, 70.0, 70.0},       // stream_width in km
    {5.0, 5.0, 5.0},           // upstream till phi
    {5.0, 3.0, 8.0}            // downstream till phi
  }

};


IcePSTexModel::IcePSTexModel(IceGrid &g, IceType *i)
  : IceEISModel(g,i) {  // do almost nothing;
                        // derived classes must have constructors
  expername = 'A';  // starts with something close to this, anyway
}


PetscErrorCode IcePSTexModel::initFromOptions() {
  PetscErrorCode      ierr;

  exper_chosen = -1;
  for (int j=0; j<Nexpers; j++) {
    PetscTruth  optionset;
    char optionname[20] = "-";
    strcat(optionname,exper_names[j]);
    // ierr = verbPrintf(2,grid.com,"optionname = %s\n",optionname);  CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, optionname, &optionset);
      CHKERRQ(ierr);
    if (optionset == PETSC_TRUE) {
      if (exper_chosen >= 0) {
        SETERRQ(1,"Only one experiment name option allowed for IcePSTexModel.");
      } else {
        exper_chosen = j;
        strcpy(exper_chosen_name, exper_names[j]);
      }
    }
  }
  if (exper_chosen < 0)
    SETERRQ(2,"Unrecognized experiment name for IcePSTexModel.\n"
              "  An experiment name option like '-P2' must be chosen.");

  if (exper_chosen <= 1) {
    useSSAVelocity = PETSC_FALSE;
    doSuperpose = PETSC_FALSE;
    doPlasticTill = PETSC_FALSE;
  } else {
    useSSAVelocity = PETSC_TRUE;
    doSuperpose = PETSC_TRUE;
    doPlasticTill = PETSC_TRUE;
  }
  
  // these are different from EISMINT II conventions
  updateHmelt = PETSC_TRUE;
  includeBMRinContinuity = PETSC_TRUE;

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);  
  
  ierr = verbPrintf(2,grid.com, 
    "setting up PST (Plastic till Stream w Thermocoupling) experiment '%s' ...\n",
    exper_chosen_name); CHKERRQ(ierr);

  ierr = setBedElev(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
         "bed topography stored by IcePSTexModel ...\n");
         CHKERRQ(ierr);

  ierr = setTillPhi(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
         "map of phi = (till friction angle) stored by IcePSTexModel ...\n");
         CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, 
     "running PST (Plastic till Stream w Thermocoupling) experiment '%s' ...\n",
    exper_chosen_name); CHKERRQ(ierr);
  return 0;
}


int IcePSTexModel::sectorNumber(const PetscScalar x, const PetscScalar y) {
  if (x > 0.0) {
    if (y < x)
      return 0;
    else
      return 1;
  } else {
    if (y < 0.0)
      return 2;
    else 
      return 1;
  }
}


PetscInt IcePSTexModel::setBedElev() {
  PetscErrorCode ierr;
  PetscScalar **b;
  
  // half-width of EISMINT II computational domain
  const PetscScalar    L = 750.0e3,
                       offset = 100.0e3,
                       width = 200.0e3,  // trough width = 200km
                       plateau = 2000.0;
  const PetscScalar    dx = grid.dx, dy = grid.dy;

  ierr = VecSet(vbed, plateau); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt m=0; m<3; m++) { // go through sectors
    PetscScalar    drop = exper_params[exper_chosen][0][m];
    PetscScalar    slope = drop / (L - offset);
    //ierr = verbPrintf(2,grid.com,
    //         "setting bed for sector %d with drop %f and angle %f\n",
    //         m,drop,stream_angle[m]); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar x = -grid.Ly + dy * j,  // note reversal
                          y = -grid.Lx + dx * i;
        if (m == sectorNumber(x,y)) {
          const PetscScalar sinrot = sin((pi/180.0)*stream_angle[m]),
                            cosrot = cos((pi/180.0)*stream_angle[m]);
          const PetscScalar x_rot =  cosrot * x + sinrot * y,
                            y_rot = -sinrot * x + cosrot * y;
          if ( (x_rot > offset) && (y_rot < width / 2.0) && (y_rot > -width / 2.0) ) {
            b[i][j] = plateau - slope * (x_rot - offset) * cos(pi * y_rot / width);
          }
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  return 0;
}




PetscInt IcePSTexModel::setTillPhi() {
  PetscErrorCode ierr;
  PetscScalar **phi;
  
  const PetscScalar    offset = 100.0e3,
                       stream_change = 400.0e3;
  const PetscScalar    dx = grid.dx, dy = grid.dy;

  ierr = VecSet(vtillphi, DEFAULT_PHI_STRONG); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vtillphi, &phi); CHKERRQ(ierr);
  for (PetscInt m=0; m<3; m++) { // go through sectors
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar x = -grid.Ly + dy * j,  // note reversal
                          y = -grid.Lx + dx * i;
        if (m == sectorNumber(x,y)) {
          const PetscScalar width = exper_params[exper_chosen][1][m] * 1000.0;
          const PetscScalar sinrot = sin((pi/180.0)*stream_angle[m]),
                            cosrot = cos((pi/180.0)*stream_angle[m]);
          const PetscScalar x_rot =  cosrot * x + sinrot * y,
                            y_rot = -sinrot * x + cosrot * y;
          if ( (x_rot > offset) && (y_rot < width / 2.0) && (y_rot > -width / 2.0) ) {
            if ((x_rot - offset) > stream_change) {
              phi[i][j] = exper_params[exper_chosen][3][m]; // downstream value
            } else {
              phi[i][j] = exper_params[exper_chosen][2][m]; // upstream value
            }
          }
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vtillphi, &phi); CHKERRQ(ierr);

  return 0;
}


//! Replace IceModel default summary.  This one gives velocities everywhere and in stream.
PetscErrorCode IcePSTexModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

// FIXME!!

  PetscErrorCode ierr;

  PetscScalar     **H, **ubar, **vbar, **ub, **vb;
  // these are in MKS; with "g" are global across all processors
  PetscScalar     maxcbar = 0.0, 
                  avcbar = 0.0, avcbupstream = 0.0, avcbdownstream = 0.0, 
                  Nhaveice = 0.0, Nupstream = 0.0, Ndownstream;
  PetscScalar     gmaxcbar, gavcbar, gavcbupstream, gavcbdownstream,
                  gNhaveice, gNupstream, gNdownstream;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvb, &vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        Nhaveice += 1.0;
        const PetscScalar cbar = sqrt( PetscSqr(ubar[i][j])
                                       + PetscSqr(vbar[i][j]) );
        if (cbar > maxcbar)  maxcbar = cbar;
        avcbar += cbar;
/*
        const PetscInt code = tillRegionCode(i,j);
        if ((code == 2) || (code == 3)) {
          const PetscScalar cb = sqrt( PetscSqr(ub[i][j])+ PetscSqr(vb[i][j]));
          if (code == 2) {
            Nupstream += 1.0;
            avcbupstream += cb;
          } else {
            Ndownstream += 1.0;
            avcbdownstream += cb;
          }
        }
*/
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vub, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvb, &vb); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxcbar, &gmaxcbar, grid.com); CHKERRQ(ierr);
  
  ierr = PetscGlobalSum(&Nhaveice, &gNhaveice, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&Nupstream, &gNupstream, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&Ndownstream, &gNdownstream, grid.com); CHKERRQ(ierr);

  ierr = PetscGlobalSum(&avcbar, &gavcbar, grid.com); CHKERRQ(ierr);
  if (gNhaveice > 0)   gavcbar = gavcbar / gNhaveice;
  else                 gavcbar = 0.0;  // degenerate case
    
  ierr = PetscGlobalSum(&avcbupstream, &gavcbupstream, grid.com); CHKERRQ(ierr);
  if (gNupstream > 0)  gavcbupstream = gavcbupstream / gNupstream;
  else                 gavcbupstream = 0.0;  // degenerate case

  ierr = PetscGlobalSum(&avcbdownstream, &gavcbdownstream, grid.com); CHKERRQ(ierr);
  if (gNdownstream > 0)  gavcbdownstream = gavcbdownstream / gNdownstream;
  else                   gavcbdownstream = 0.0;  // degenerate case

  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "P         YEAR:     ivol   iarea   meltf   maxcbar    avcbar   avcbUpS avcbDownS\n");
    ierr = verbPrintf(2,grid.com,
      "U        years 10^6_km^3 10^6_km^2 (none)      m/a       m/a       m/a       m/a\n");
  } else {
    ierr = verbPrintf(2,grid.com, 
      "S %12.5f: %8.5f %7.4f %7.4f %9.2f %9.3f %9.3f %9.3f\n",
         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, meltfrac,
         gmaxcbar*secpera, gavcbar*secpera, 
         gavcbupstream*secpera, gavcbdownstream*secpera); CHKERRQ(ierr);
  }
  return 0;
}
