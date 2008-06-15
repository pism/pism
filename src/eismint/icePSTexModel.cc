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

const PetscScalar 
  DEFAULT_PHI_STRONG = 15.0,  // till friction angle outside of stream
  stream_offset = 100.0e3,    // distance from grid center to start of stream
  stream_length = 650.0e3,    // length used in computing slope
  stream_change = 400.0e3,    // distance along stream at which  change from
                              //     'upstream' to 'downstream' till friction
                              //     angle occurs
  xi_slop = 0.15,             // fractions by which reduced till friction
  eta_slop = 0.5;            //     angle extends outside of stream; see
                              //     setTillPhi() and inStreamNbhd()

const int Nexpers = 6,
          Nstreams = 4,     // 4 streams for most experiments, 3 for P2
          NAME_LENGTH = 10;

struct ExperDescription {
  char         name[NAME_LENGTH];
  PetscScalar  bed_end_depth[4];   // drop for stream (see stream_length)
  PetscScalar  stream_width[4];    // width in km
  PetscScalar  upstream_phi[4];    // till friction angle in degrees
  PetscScalar  downstream_phi[4];
};

ExperDescription e[Nexpers] = {

  { "P0A",
    {0.0, 0.0, 0.0, 0.0}, 
    {70.0, 70.0, 70.0, 70.0}, // INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // INACTIVE
  },

  { "P0I",
    {1000.0, 1000.0, 1000.0, 1000.0},
    {70.0, 70.0, 70.0, 70.0}, // INACTIVE
    {5.0, 5.0, 5.0, 5.0},     // INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // INACTIVE
  },

  { "P1",                     // default experiment
    {0.0, 0.0, 0.0, 0.0},
    {70.0, 30.0, 100.0, 50.0},
    {5.0, 5.0, 5.0, 5.0},
    {5.0, 5.0, 5.0, 5.0}
  },

  { "P2",                     // only three streams, so last entry inactive
    {0.0, 0.0, 0.0, 0.0},     // last one INACTIVE
    {70.0, 70.0, 70.0, 70.0}, // last *three* INACTIVE; constant width
    {5.0, 5.0, 5.0, 5.0},     // last one INACTIVE
    {5.0, 5.0, 5.0, 5.0}      // last one INACTIVE
  },

  { "P3",
    {1000.0, 1000.0, 1000.0, 1000.0},
    {70.0, 30.0, 100.0, 50.0},
    {5.0, 5.0, 5.0, 5.0},
    {5.0, 5.0, 5.0, 5.0}
  },

  { "P4",
    {0.0, 0.0, 0.0, 0.0},
    {70.0, 70.0, 70.0, 70.0},
    {5.0, 5.0, 5.0, 5.0},
    {2.0, 3.0, 8.0, 10.0}
  }

};

PetscScalar stream_angle_P2[3] = {0.0, 100.0, 225.0};  // degrees


IcePSTexModel::IcePSTexModel(IceGrid &g, IceType *i)
  : IceEISModel(g,i) {  // do almost nothing; derived need constructors
  expername = 'A';      // mostly close to this, anyway
}


PetscErrorCode IcePSTexModel::initFromOptions() {
  PetscErrorCode      ierr;

  exper_chosen = -1;
  for (int j=0; j<Nexpers; j++) {
    PetscTruth  optionset;
    char optionname[20] = "-";
    strcat(optionname,e[j].name);
    ierr = PetscOptionsHasName(PETSC_NULL, optionname, &optionset);
      CHKERRQ(ierr);
    if (optionset == PETSC_TRUE) {
      if (exper_chosen >= 0) {
        SETERRQ(1,"Only one experiment name option allowed for IcePSTexModel.");
      } else {
        exper_chosen = j;
        strcpy(exper_chosen_name,e[j].name);
      }
    }
  }
  if (exper_chosen < 0)
    SETERRQ(2,"Unrecognized experiment name for IcePSTexModel.\n"
              "  An experiment name option like '-P2' must be chosen.");

  ierr = verbPrintf(2,grid.com, 
    "setting up PST (Plastic till Stream w Thermocoupling) experiment '%s' ...\n",
    exper_chosen_name); CHKERRQ(ierr);

  ierr = IceEISModel::initFromOptions(); CHKERRQ(ierr);  

  // different from EISMINT II conventions (even for P0A and P0I)
  updateHmelt = PETSC_TRUE;
  includeBMRinContinuity = PETSC_TRUE;
  transformForSurfaceGradient = PETSC_TRUE;

  if (exper_chosen <= 1) { // P0A and P0I are nonsliding SIA
    useSSAVelocity = PETSC_FALSE;
    doSuperpose = PETSC_FALSE;
    doPlasticTill = PETSC_FALSE;
  } else {
    // these options equiv to "-ssa -super -plastic"
    useSSAVelocity = PETSC_TRUE;
    doSuperpose = PETSC_TRUE;
    doPlasticTill = PETSC_TRUE;
    useConstantTillPhi = PETSC_FALSE;
  }  

  ierr = setBedElev(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, "bed topography stored ... "); CHKERRQ(ierr);

  ierr = setTillPhi(); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
         "map of phi = (till friction angle) stored\n"); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com, 
     "running PST experiment '%s' ...\n",
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


//! Say whether we are in the stream (strictly or not), and give local coords.
bool IcePSTexModel::inStreamNbhd(bool strictly_in_stream,
                            const PetscScalar angle, const PetscScalar width,
                            const PetscScalar x, const PetscScalar y,
                            PetscScalar &x_loc, PetscScalar &y_loc) {
  const PetscScalar sinrot = sin(angle),
                    cosrot = cos(angle);
  x_loc =  cosrot * x + sinrot * y - stream_offset;
  y_loc = -sinrot * x + cosrot * y;
  if (strictly_in_stream) {
    if ( (x_loc > 0.0) && (fabs(y_loc) < width / 2.0) )
      return true;
    else
      return false;
  } else {
    if ( (x_loc > - xi_slop * stream_change)
         && (fabs(y_loc) < (1.0 + eta_slop) * (width / 2.0)) )
      return true;
    else
      return false;
  }
}


//! Say whether we are strictly in the stream, and give local coords.
bool IcePSTexModel::inStream(const PetscScalar angle, const PetscScalar width,
                            const PetscScalar x, const PetscScalar y,
                            PetscScalar &x_loc, PetscScalar &y_loc) {
  return inStreamNbhd(true, angle,width, x,y, x_loc,y_loc);
}


PetscErrorCode IcePSTexModel::setBedElev() {
  PetscErrorCode ierr;
  PetscScalar **b;
  
  const PetscScalar    width = 200.0e3,  // trough width = 200km; not the same
                                         //   as stream width
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
        PetscScalar drop = e[exper_chosen].bed_end_depth[m],
                    slope = drop / stream_length;
        if (inStream((pi/2.0)*m,width,x,y,x_loc,y_loc))
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
  if (PetscAbs(eta) <= 1.0 - eta_slop)
    lambda = 1.0;
  else if (PetscAbs(eta) < 1.0 + eta_slop)
    lambda = 0.5 - 0.5 * sin((pi/2.0) * (PetscAbs(eta) - 1.0) / eta_slop);

  if (x > stream_change)
    return DOWN * lambda + STRONG * (1.0 - lambda); // downstream value
  else { // f(xi) is for upstream part only
    PetscScalar f = STRONG;
    if (xi >= xi_slop)
      f = UP;
    else if (xi > - xi_slop) {
      const PetscScalar fav = 0.5 * (STRONG + UP);
      f = fav - 0.5 * (STRONG - UP) * sin((pi/2.0) * (xi / xi_slop));
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
    return 0;  // nothing further for P0A and P0I

  ierr = DAVecGetArray(grid.da2, vtillphi, &phi); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar x = -grid.Ly + dy * j,  // note reversal
                        y = -grid.Lx + dx * i;
      if (exper_chosen == 3) { // experiment P2
        const PetscScalar width = e[exper_chosen].stream_width[0] * 1000.0;
        for (PetscInt m=0; m<3; m++) {
          if (inStreamNbhd(false,(pi/180.0) * stream_angle_P2[m],width,x,y,x_loc,y_loc))
            phi[i][j] = phiLocal(width,x_loc,y_loc,DEFAULT_PHI_STRONG,
                                 e[exper_chosen].upstream_phi[m],
                                 e[exper_chosen].downstream_phi[m]);
        }
      } else {
        for (PetscInt m=0; m<4; m++) { // go through four sectors
          const PetscScalar width = e[exper_chosen].stream_width[m] * 1000.0;
          if (inStreamNbhd(false,(pi/2.0)*m,width,x,y,x_loc,y_loc)) {
            phi[i][j] = phiLocal(width,x_loc,y_loc,DEFAULT_PHI_STRONG,
                                 e[exper_chosen].upstream_phi[m],
                                 e[exper_chosen].downstream_phi[m]);
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
          const PetscScalar width = e[exper_chosen].stream_width[0] * 1000.0;
          for (int m=0; m<3; m++) {
            if (inStream((pi/180.0) * stream_angle_P2[m],width,x,y,x_loc,y_loc)) {
              if (r > stream_change) {
                areadown[m] += darea;
                avcdown[m] += cbar * darea;
              } else {
                areaup[m] += darea;
                avcup[m] += cbar * darea;
              }
            }
          }
        } else {
          for (int m=0; m<4; m++) {
            const PetscScalar width = e[exper_chosen].stream_width[m] * 1000.0;
            if (inStream((pi/2.0)*m,width,x,y,x_loc,y_loc)) {
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
      "++++++++++++++++++++++++++ WIDEN SCREEN TO AT LEAST %d CHARS FOR IcePSTexModel OUTPUT"
      " +++++++++++++++++++++++++++\n",
      (exper_chosen == 3) ? 96 : 114); CHKERRQ(ierr);
    if (exper_chosen == 3) {
      ierr = verbPrintf(2,grid.com,
        "P         YEAR:     ivol   iarea   maxcbar"
        "    avup0   avdwn0    avup1   avdwn1    avup2   avdwn2\n");
      ierr = verbPrintf(2,grid.com,
        "U        years 10^6_km^3 10^6_km^2     m/a"
        "      m/a      m/a      m/a      m/a      m/a      m/a\n");
    } else {
      ierr = verbPrintf(2,grid.com,
        "P         YEAR:     ivol   iarea   maxcbar"
        "    avup0   avdwn0    avup1   avdwn1    avup2   avdwn2    avup3   avdwn3\n");
      ierr = verbPrintf(2,grid.com,
        "U        years 10^6_km^3 10^6_km^2     m/a"
        "      m/a      m/a      m/a      m/a      m/a      m/a      m/a      m/a\n");
    }
  } else { // run time summary line
    if (exper_chosen == 3) {
      ierr = verbPrintf(2,grid.com, 
        "S %12.5f: %8.5f %7.4f %9.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6,
         gmaxcbarALL*secpera,
         gavcup[0]*secpera, gavcdown[0]*secpera,
         gavcup[1]*secpera, gavcdown[1]*secpera,
         gavcup[2]*secpera, gavcdown[2]*secpera); CHKERRQ(ierr);
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
  }

  return 0;
}

