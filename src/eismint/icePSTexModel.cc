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

PetscScalar stream_angle_P2[3] = {0.0, 100.0, 225.0};  // degrees

const PetscScalar DEFAULT_PHI_STRONG = 20.0;

const PetscScalar stream_length = 650.0e3;

const PetscScalar stream_change = 400.0e3;

PetscScalar exper_params[Nexpers][4][4] = {

  {// P0A
    {0.0, 0.0, 0.0, 0.0},     // bed depth at end
    {70.0, 70.0, 70.0, 70.0}, // stream_width in km; INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi; INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // downstream till phi; INACTIVE
  },

  {// P0I
    {1000.0, 1000.0, 1000.0, 1000.0},    // bed depth at end
    {70.0, 70.0, 70.0, 70.0}, // stream_width in km; INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi; INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // downstream till phi; INACTIVE
  },

  {// P1: default
    {0.0, 0.0, 0.0, 0.0},     // bed depth at end
    {70.0, 30.0, 100.0, 50.0}, // stream_width in km
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi
    {5.0, 5.0, 5.0, 5.0}      // downstream till phi
  },

  {// P2: only three streams, so last entry inactive
    {0.0, 0.0, 0.0, 0.0},     // bed depth at end
    {70.0, 70.0, 70.0, 70.0}, // stream_width in km; INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi
    {5.0, 5.0, 5.0, 5.0}      // downstream till phi
  },

  {// P3
    {1000.0, 1000.0, 1000.0, 1000.0},    // bed depth at end
    {70.0, 30.0, 100.0, 50.0}, // stream_width in km
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi
    {5.0, 5.0, 5.0, 5.0}      // downstream till phi
  },

  {// P4
    {0.0, 0.0, 0.0, 0.0},     // bed depth at end
    {70.0, 70.0, 70.0, 70.0}, // stream_width in km
    {5.0, 5.0, 5.0, 5.0},     // upstream till phi
    {2.0, 3.0, 8.0, 10.0}     // downstream till phi
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

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);  

  if (exper_chosen <= 1) {
    useSSAVelocity = PETSC_FALSE;
    doSuperpose = PETSC_FALSE;
    doPlasticTill = PETSC_FALSE;
  } else {
    // these options mean "-ssa -super -plastic" already
    useSSAVelocity = PETSC_TRUE;
    doSuperpose = PETSC_TRUE;
    doPlasticTill = PETSC_TRUE;
    useConstantTillPhi = PETSC_FALSE;
  }
  
  // these are different from EISMINT II conventions
  updateHmelt = PETSC_TRUE;
  includeBMRinContinuity = PETSC_TRUE;
  
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


int IcePSTexModel::sectorNumberP2(const PetscScalar x, const PetscScalar y) {
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


//! Echo four-way stream number if in stream, or -1 if not in stream, and give local coords in stream.
int IcePSTexModel::inStream(const int m, const PetscScalar width,
                            const PetscScalar x, const PetscScalar y,
                            PetscScalar &x_loc, PetscScalar &y_loc) {
  const PetscScalar offset = 100.0e3,
                    sinrot = sin(0.5 * pi * m),
                    cosrot = cos(0.5 * pi * m);
  x_loc =  cosrot * x + sinrot * y - offset;
  y_loc = -sinrot * x + cosrot * y;
  if ( (x_loc > 0.0) && (y_loc < width / 2.0) && (y_loc > -width / 2.0) )
    return m;
  else
    return -1;
}


//! Return three-way stream number (for experiment P2) or -1 if not in stream.
int IcePSTexModel::inStreamP2(const PetscScalar width,
                              const PetscScalar x, const PetscScalar y,
                              PetscScalar &x_loc, PetscScalar &y_loc) {
  const int m = sectorNumberP2(x,y);
  const PetscScalar offset = 100.0e3,
                    sinrot = sin((pi/180.0) * stream_angle_P2[m]),
                    cosrot = cos((pi/180.0) * stream_angle_P2[m]);
  x_loc =  cosrot * x + sinrot * y - offset;
  y_loc = -sinrot * x + cosrot * y;
  if ( (x_loc > 0.0) && (y_loc < width / 2.0) && (y_loc > -width / 2.0) )
    return m;
  else
    return -1;
}


PetscErrorCode IcePSTexModel::setBedElev() {
  PetscErrorCode ierr;
  PetscScalar **b;
  
  const PetscScalar    width = 200.0e3,  // trough width = 200km
                       plateau = 2000.0;
  const PetscScalar    dx = grid.dx, dy = grid.dy;
  PetscScalar x_loc, y_loc;

  ierr = VecSet(vbed, plateau); CHKERRQ(ierr);

  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar x = -grid.Ly + dy * j,  // note reversal
                        y = -grid.Lx + dx * i;
      // note we treat exper P2 like others; it is flat anyway (slope=0)
      for (PetscInt m=0; m<4; m++) {
        PetscScalar drop = exper_params[exper_chosen][0][m],
                    slope = drop / stream_length;
        if (m == inStream(m,width,x,y,x_loc,y_loc))
          b[i][j] = plateau - slope * x_loc * cos(pi * y_loc / width);
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  return 0;
}


//! Compute the till friction angle in local strip coordinates.
PetscScalar IcePSTexModel::phiLocal(const PetscScalar width,
              const PetscScalar x, const PetscScalar y,
              const PetscScalar STRONG, 
              const PetscScalar UP, const PetscScalar DOWN) {

  const PetscScalar eta   = y / (width/2.0),   // normalized local y
                    xi    = x / stream_change; // normalized local x

  // compute lambda(eta) which is even and in [0,1]
  PetscScalar lambda = 0.0; // for big eta
  if (PetscAbs(eta) <= 0.5)
    lambda = 1.0;
  else if (PetscAbs(eta) < 1.5)
    lambda = 0.5 - 0.5 * sin(pi * (PetscAbs(eta) - 1.0));

  if (x > stream_change)
    return DOWN * lambda + STRONG * (1.0 - lambda); // downstream value
  else { // f(xi) is for upstream part only
    PetscScalar f = STRONG;
    if (xi >= 0.25)
      f = UP;
    else if (xi > -0.25) {
      const PetscScalar fav = 0.5 * (STRONG + UP);
      f = fav - 0.5 * sin(2.0 * pi * xi) * (STRONG - UP);
    }
    return f * lambda + STRONG * (1.0 - lambda); // upstream value
  }
}



PetscErrorCode IcePSTexModel::setTillPhi() {
  PetscErrorCode ierr;
  PetscScalar **phi;
  
  const PetscScalar    dx = grid.dx, dy = grid.dy;
  PetscScalar          x_loc, y_loc;

  ierr = VecSet(vtillphi, DEFAULT_PHI_STRONG); CHKERRQ(ierr);

  if (exper_chosen <= 1)
    return 0;  // done for P0A and P0I

  ierr = DAVecGetArray(grid.da2, vtillphi, &phi); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar x = -grid.Ly + dy * j,  // note reversal
                        y = -grid.Lx + dx * i;
      // FIXME:  need to return smoothed phi even if not *IN* stream
      //   but merely near the stream; create nearStream, nearStreamP2
      if (exper_chosen == 3) { // experiment P2
        const PetscScalar width = exper_params[exper_chosen][1][0] * 1000.0;
        int m = inStreamP2(width,x,y,x_loc,y_loc);
        if (m >= 0)
          phi[i][j] = phiLocal(width,x_loc,y_loc,DEFAULT_PHI_STRONG,
                               exper_params[exper_chosen][2][m],
                               exper_params[exper_chosen][3][m]);
      } else {
        for (PetscInt m=0; m<4; m++) { // go through four sectors
          const PetscScalar width = exper_params[exper_chosen][1][m] * 1000.0;
          if (m == inStream(m,width,x,y,x_loc,y_loc)) {
            phi[i][j] = phiLocal(width,x_loc,y_loc,DEFAULT_PHI_STRONG,
                                 exper_params[exper_chosen][2][m],
                                 exper_params[exper_chosen][3][m]);
          }          
        }
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vtillphi, &phi); CHKERRQ(ierr);

  return 0;
}


//! Replace IceModel default summary.  This one gives velocities in streams.
PetscErrorCode IcePSTexModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

  PetscErrorCode ierr;
  PetscScalar     **H, **ubar, **vbar;
  // these are in MKS units; with "g" are global across all processors
  PetscScalar     maxcbarALL = 0.0, 
                  areaup[4] = {0.0, 0.0, 0.0, 0.0},
                  avcup[4] = {0.0, 0.0, 0.0, 0.0},
                  areadown[4] = {0.0, 0.0, 0.0, 0.0},
                  avcdown[4] = {0.0, 0.0, 0.0, 0.0};
  PetscScalar     gmaxcbarALL, gareaup[4], gavcup[4], gareadown[4], gavcdown[4];
  PetscScalar     x_loc, y_loc;
  const PetscScalar darea = grid.dx * grid.dy;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0) {
        const PetscScalar cbar = sqrt( PetscSqr(ubar[i][j])
                                       + PetscSqr(vbar[i][j]) );
        const PetscScalar x = -grid.Ly + grid.dy * j,  // note reversal
                          y = -grid.Lx + grid.dx * i,
                          r = sqrt(PetscSqr(x) + PetscSqr(y));
        if (cbar > maxcbarALL)  maxcbarALL = cbar;
        if (exper_chosen == 3) { // exper P2
          const PetscScalar width = exper_params[exper_chosen][1][0] * 1000.0;
          int m = inStreamP2(width,x,y,x_loc,y_loc);
          if (m >= 0) {
            if (r > stream_change) {
              areadown[m] += darea;
              avcdown[m] += cbar * darea;
            } else {
              areaup[m] += darea;
              avcup[m] += cbar * darea;
            }
          }
        } else {
          for (int m=0; m<4; m++) {
            const PetscScalar width = exper_params[exper_chosen][1][m] * 1000.0;
            if (m == inStream(m,width,x,y,x_loc,y_loc)) {
              if (x_loc > stream_change) {
                areadown[m] += darea;
                avcdown[m] += cbar * darea;
              } else {
                areaup[m] += darea;
                avcup[m] += cbar * darea;
              }
            }
          }
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubar, &ubar); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);

  // globalize and actually compute averages
  ierr = PetscGlobalMax(&maxcbarALL, &gmaxcbarALL, grid.com); CHKERRQ(ierr);
  for (PetscInt m=0; m<4; m++) {
    ierr = PetscGlobalSum(&areaup[m], &gareaup[m], grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&avcup[m], &gavcup[m], grid.com); CHKERRQ(ierr);
    if (gareaup[m] > 0.0) {
      gavcup[m] = gavcup[m] / gareaup[m];
    } else {
      gavcup[m] = 0.0;
    }
    ierr = PetscGlobalSum(&areadown[m], &gareadown[m], grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&avcdown[m], &gavcdown[m], grid.com); CHKERRQ(ierr);
    if (gareadown[m] > 0.0) {
      gavcdown[m] = gavcdown[m] / gareadown[m];
    } else {
      gavcdown[m] = 0.0;
    }
  }

  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "++++++++++++++++++++++++++ WIDEN SCREEN TO AT LEAST 114 CHARS FOR IcePSTexModel OUTPUT"
      " +++++++++++++++++++++++++++\n"); CHKERRQ(ierr);
    if (exper_chosen == 3) {
      ierr = verbPrintf(2,grid.com, 
        "\nWARNING: columns 'avup3' and 'avdwn3' should be ignored for experiment P2.\n\n");
        CHKERRQ(ierr);
    }
    ierr = verbPrintf(2,grid.com,
      "P         YEAR:     ivol   iarea   maxcbar"
      "    avup0   avdwn0    avup1   avdwn1    avup2   avdwn2    avup3   avdwn3\n");
    ierr = verbPrintf(2,grid.com,
      "U        years 10^6_km^3 10^6_km^2     m/a"
      "      m/a      m/a      m/a      m/a      m/a      m/a      m/a      m/a\n");
  } else {
    ierr = verbPrintf(2,grid.com, 
      "S %12.5f: %8.5f %7.4f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6,
         gmaxcbarALL*secpera,
         gavcup[0]*secpera, gavcdown[0]*secpera,
         gavcup[1]*secpera, gavcdown[1]*secpera,
         gavcup[2]*secpera, gavcdown[2]*secpera,
         gavcup[3]*secpera, gavcdown[3]*secpera); CHKERRQ(ierr);
  }

  return 0;
}

