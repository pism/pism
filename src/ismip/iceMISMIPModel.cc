// Copyright (C) 2008 Ed Bueler
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
#include <petscda.h>
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "iceMISMIPModel.hh"


/* 
typical run will be
  $ pisms -mismip 2b -initials EBU -mode 1 -run 3
will automatically run until steady state criterion or 3e4 years, and save NetCDF file
  EBU1_2b_M1_A3.nc
and text files
  EBU1_2b_M1_A3_t
  EBU1_2b_M1_A3_ss
  EBU1_2b_M1_A3_f

alternate run might be
  $ pisms -mismip 2b -initials EBU -mode 1 -run 3 -super
will automatically run until steady state criterion or 3e4 years, and save NetCDF file
  EBU2_2b_M1_A3.nc
and text files
  EBU2_2b_M1_A3_t
  EBU2_2b_M1_A3_ss
  EBU2_2b_M1_A3_f
*/


PetscScalar MISMIPIce::flow(const PetscScalar stress, const PetscScalar temp,
                            const PetscScalar pressure) const {
  return A_MISMIP * pow(stress,n-1);
}


PetscScalar MISMIPIce::effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const  {
  const PetscScalar nn = (PetscScalar) n;
  const PetscScalar alpha = 0.5 * PetscSqr(u_x) + 0.5 * PetscSqr(v_y)
                             + 0.5 * PetscSqr(u_x + v_y) + 0.25 * PetscSqr(u_y + v_x);
  return 0.5 * B_MISMIP * pow(regularization + alpha, - (nn - 1.0) / (2.0 * nn));
}


PetscScalar MISMIPIce::effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const  {
  // DESPITE NAME, does *not* return effective viscosity; returns viscosity times thickness.
  // NOTE: temp and pressure args to effectiveViscosity() ignored
  return H * effectiveViscosity(regularization, u_x, u_y, v_x, v_y, 0.0, 0.0);
}


MISMIPIce::MISMIPIce() {
  rho = 900.0;
  n = 3;
  setA(4.6416e-24);
}


PetscErrorCode MISMIPIce::setA(const PetscScalar myA) {
  A_MISMIP = myA;
  B_MISMIP = pow(A_MISMIP, - 1.0 / ((PetscScalar) n));
  return 0;
}


PetscScalar MISMIPIce::getA() {
  return A_MISMIP;
}


PetscErrorCode MISMIPIce::printInfo(const int thresh, MPI_Comm com) {
  PetscErrorCode ierr;
  ierr = verbPrintf(thresh, com, 
    "Using MISMIP ice with rho = %5.4f, n = %d, and A=%5.4e.  (Note grav = %5.4f.)\n",
    rho, n, A_MISMIP, grav); CHKERRQ(ierr);
  return 0;
}


MISMIPBasalType::MISMIPBasalType(const PetscScalar m, const PetscScalar C,
                                 const PetscScalar regularize) {
  m_MISMIP = m;
  C_MISMIP = C;
  regularize_MISMIP = regularize;
}


PetscErrorCode MISMIPBasalType::printInfo(const int thresh, MPI_Comm com) {
  PetscErrorCode ierr;
  if (m_MISMIP == 1.0) {
    ierr = verbPrintf(thresh, com, 
      "Using  tau_b = - C u  sliding model for MISMIP, with C=%5.4e.\n", C_MISMIP); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(thresh, com, 
      "Using  tau_b = - C (|u|^2 + eps^2)^{(m-1)/2} u  model for MISMIP, with m=%5.4f, C=%5.4e,\n"
      "  and eps = %5.4f m/a.\n",
      m_MISMIP, C_MISMIP, regularize_MISMIP * secpera); CHKERRQ(ierr);
  }
  return 0;
}


PetscScalar MISMIPBasalType::drag(PetscScalar coeff, PetscScalar tauc,
                                   PetscScalar vx, PetscScalar vy) {
  if (m_MISMIP == 1.0) {
    return C_MISMIP;
  } else {
    const PetscScalar magsliding = PetscSqr(vx) + PetscSqr(vy) + PetscSqr(regularize_MISMIP);
    return C_MISMIP * pow(magsliding, (m_MISMIP - 1.0) / 2.0);
  }
}



IceMISMIPModel::IceMISMIPModel(IceGrid &g, IceType &i, MISMIPIce &mismip_i) : 
      IceModel(g, i), mismip_ice(mismip_i) {

  // this just fills the data members with non-junk; see setFromOptions() for decent values
  exper = 1;
  sliding = 'a';
  gridmode = 1;
  runindex = 1;
  runtimeyears = 3.0e4;
  strcpy(initials,"ABC");
}


PetscErrorCode IceMISMIPModel::setFromOptions() {
  PetscErrorCode ierr;

  // from Table 4
  const PetscScalar Aexper1or2[10] = {0.0, // zero position not used
                        4.6416e-24,  2.1544e-24,  1.0e-24,
                        4.6416e-25,  2.1544e-25,  1.0e-25,
                        4.6416e-26,  2.1544e-26,  1.0e-26};
  // from Table 5
  const PetscScalar timeexper3a[14] = {0.0, // zero position not used
                        3.0e4, 1.5e4, 1.5e4, 
                        1.5e4, 1.5e4, 3.0e4,
                        3.0e4, 1.5e4, 1.5e4,
                        3.0e4, 3.0e4, 3.0e4,
                        1.5e4};
  const PetscScalar Aexper3a[14] = {0.0, // zero position not used
                        3.0e-25, 2.5e-25, 2.0e-25,
                        1.5e-25, 1.0e-25, 5.0e-26,
                        2.5e-26, 5.0e-26, 1.0e-25,
                        1.5e-25, 2.0e-25, 2.5e-25,
                        3.0e-25};
  // from Table 6
  const PetscScalar timeexper3b[16] = {0.0, // zero position not used
                        3.0e4, 1.5e4, 1.5e4,
                        1.5e4, 1.5e4, 1.5e4,
                        1.5e4, 3.0e4, 1.5e4,
                        1.5e4, 1.5e4, 1.5e4,
                        1.5e4, 3.0e4, 1.5e4};   //  15th VALUE LABELED AS 16 IN Table 6 !?
  const PetscScalar Aexper3b[16] = {0.0, // zero position not used
                        1.6e-24, 1.4e-24, 1.2e-24,
                        1.0e-24, 8.0e-25, 6.0e-25,
                        4.0e-25, 2.0e-25, 4.0e-25,
                        6.0e-25, 8.0e-25, 1.0e-24,
                        1.2e-24, 1.4e-24, 1.6e-24};   //  15th VALUE LABELED AS 16 IN Table 6 !?

  // read option    -mismip [1|1a|1b|2|2a|2b|3|3a|3b]
  char Ee[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-mismip", Ee, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (strlen(Ee) == 0) {
    exper = 1;
    sliding = 'a';
  } else if (strlen(Ee) == 1) {
    if ((Ee[0] < '1') || (Ee[0] > '3')) {
      SETERRQ(1,"IceMISMIPModel ERROR:  first character of string Ee in '-mismip Ee' must be 1, 2, or 3");
    }
    exper = (int) Ee[0] - (int) '0';
    sliding = 'a';
  } else if (strlen(Ee) == 2) {
    if ((Ee[0] < '1') || (Ee[0] > '3')) {
      SETERRQ(1,"IceMISMIPModel ERROR:  first character of string Ee in '-mismip Ee' must be 1, 2, or 3");
    }
    exper = (int) Ee[0] - (int) '0';
    if ((Ee[1] == 'a') || (Ee[1] == 'b')) {
      sliding = Ee[1];
    } else {
      SETERRQ(2,"IceMISMIPModel ERROR:  second character of string Ee in '-mismip Ee' must be a or b");
    }
  } else {
    SETERRQ(3,"IceMISMIPModel ERROR:  string Ee in '-mismip Ee' must be of length 0, 1, or 2");
  }

  // read option    -mode [1|2|3]
  ierr = PetscOptionsGetInt(PETSC_NULL, "-mode", &gridmode, PETSC_NULL); CHKERRQ(ierr);
  if ((gridmode < 1) || (gridmode > 3)) {
    SETERRQ(8,"IceMISMIPModel ERROR:  grid mode N in '-mode N' must be 1, 2, or 3");
  }

  // read option    -run [1|..|15]
  ierr = PetscOptionsGetInt(PETSC_NULL, "-run", &runindex, PETSC_NULL); CHKERRQ(ierr);
  if (runindex < 1) {
    SETERRQ(4,"IceMISMIPModel ERROR:  run index N in '-run N' must be at least 1");
  }
  if ((exper == 1) || (exper == 2)) {
    if (runindex > 9) {
      SETERRQ(5,"IceMISMIPModel ERROR:  run index N in '-run N' must be <= 9 in experiments 1 or 2");
    }
    runtimeyears = 3.0e4;
    ierr = mismip_ice.setA(Aexper1or2[runindex]); CHKERRQ(ierr);  
  } else if (exper == 3) {
    if (sliding == 'a') {
      if (runindex > 13) {
        SETERRQ(6,"IceMISMIPModel ERROR:  run index N in '-run N' must be <= 13 in experiment 3a");
      }
      runtimeyears = timeexper3a[runindex];
      ierr = mismip_ice.setA(Aexper3a[runindex]); CHKERRQ(ierr);  
    } else if (sliding == 'b') {
      if (runindex > 15) {
        SETERRQ(7,"IceMISMIPModel ERROR:  run index N in '-run N' must be <= 15 in experiment 3b");
      }
      runtimeyears = timeexper3b[runindex];
      ierr = mismip_ice.setA(Aexper3b[runindex]); CHKERRQ(ierr);  
    } else {
      SETERRQ(99, "how did I get here?");
    }
  }

  // read option  -initials ABC
  ierr = PetscOptionsGetString(PETSC_NULL, "-initials", initials, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (strlen(initials) != 3) {
    ierr = verbPrintf(1,grid.com,"IceMISMIPModel WARNING:  Initials string should be three chars long.");
       CHKERRQ(ierr);
  }

  doTemp               = PETSC_FALSE;
  doPlasticTill        = PETSC_FALSE;
  doBedDef             = PETSC_FALSE;
  doGrainSize          = PETSC_FALSE;
  
  computeSurfGradInwardSSA = PETSC_TRUE;
  computeSIAVelocities = PETSC_FALSE;
  useSSAVelocity       = PETSC_TRUE;
  doSuperpose          = PETSC_FALSE; // FIXME: use this for model results 2?

// use value of gridmode to set grid??  set grid.Mx, grid.My from options or lack thereof

// use -y option, if given, to overwrite runtimeyears

// check on -super option;

// create name string for final file?  "EBU1_1a_M..."

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  
  return 0;
}


PetscErrorCode IceMISMIPModel::initFromOptions() {
  PetscErrorCode ierr;
  bool infileused;
  PetscTruth inFileSet, bootFileSet;

  // check if input file was used
  ierr = PetscOptionsHasName(PETSC_NULL, "-if", &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-bif", &bootFileSet); CHKERRQ(ierr);
  infileused = ((inFileSet == PETSC_TRUE) || (bootFileSet == PETSC_TRUE));

  if (infileused) {
    ierr = verbPrintf(1,grid.com, 
       "IceMISMIPModel WARNING: -if or -bif option used; not using MISMIP formulae to initialize\n");
       CHKERRQ(ierr);  
  } else { // usual case: initialize from MISMIP formulas
    ierr = verbPrintf(1,grid.com, 
              "initializing MISMIP experiment %d%c\n"
              "    [grid mode %d; run %d (A = %5.4e, run length = %5.2f years)] ... \n", 
              exper,sliding,gridmode,runindex,mismip_ice.getA(),runtimeyears); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);

    // see Table 3
    const PetscScalar REGULARIZING_VELOCITY = 0.01 / secpera;
    if (sliding == 'a') {
      basal = new MISMIPBasalType(1.0/3.0, 7.624e6, REGULARIZING_VELOCITY);
    } else if (sliding == 'b') {
      basal = new MISMIPBasalType(1.0, 7.2082e10, REGULARIZING_VELOCITY);
    } else {
      SETERRQ(99, "how did I get here?");
    }
    createBasal_done = PETSC_TRUE;

    // FIXME:  use value of gridmode to set grid??

    const PetscScalar   L = 1700.0e3;      // Horizontal half-width of grid
    // NOTE: y takes place of x!!!
    // NOTE: put calving front at 1600.0e3  ??
    ierr = grid.rescale(200.0e3, L, 4000.0); CHKERRQ(ierr);  // what to put for thickness?

    // all of these relate to models which should be turned off ...
    ierr = VecSet(vHmelt, 0.0); CHKERRQ(ierr);
    setInitialAgeYears(DEFAULT_INITIAL_AGE_YEARS);
    setConstantGrainSize(DEFAULT_GRAIN_SIZE);  // no Goldsby-Kohlstedt
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);  // no bed deformation
    ierr = VecSet(vTs, ice.meltingTemp); CHKERRQ(ierr);
    ierr = VecSet(vT, ice.meltingTemp); CHKERRQ(ierr);
    ierr = VecSet(vTb, ice.meltingTemp); CHKERRQ(ierr);

    ierr = VecSet(vAccum, 0.3/secpera); CHKERRQ(ierr);

    ierr = VecSet(vH, 10.0); CHKERRQ(ierr);
    ierr = setMISMIPBed(); CHKERRQ(ierr);

    ierr = VecSet(vMask, MASK_DRAGGING); CHKERRQ(ierr);  // note updateSurfaceElevationAndMask() will
                                                      //   set floating part to floating

    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
    
    initialized_p = PETSC_TRUE;
  }
  
  ierr = mismip_ice.printInfo(2,grid.com); CHKERRQ(ierr);
  ierr = basal->printInfo(2,grid.com); CHKERRQ(ierr);
  
  if (!isInitialized()) {
    SETERRQ(1, "ERROR: IceMISMIPModel has not been initialized!\n");
  }

  ierr = IceModel::initFromOptions(PETSC_TRUE); CHKERRQ(ierr);  // regridding can happen here

  // FIXME: report on these flags: doTemp=false, doBedDef=false, doPlasticTill=false,
  //                               doGrainSize=false,

  // check "-ssa" is set

  return 0;
}


PetscErrorCode IceMISMIPModel::setMISMIPBed() {
  PetscErrorCode ierr;
  PetscScalar          **b;

  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
      const PetscScalar jfrom0 =
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.p->My - 1)/2.0;
      const PetscScalar y = grid.p->dy * jfrom0;
      const PetscScalar xs = PetscAbs(y) / 750.0e3;  // scaled and symmetrical x coord

      if ((exper == 1) || (exper == 2)) {
        b[i][j] = 720.0 - 778.5 * xs;
      } else if (exper == 3) {
        const PetscScalar xs2 = PetscSqr(xs),
                          xs4 = PetscSqr(xs2),
                          xs6 = xs4 * xs2;
        b[i][j] = 729.0 - 2184.0 * xs2 + 1031.72 * xs4 - 151.72 * xs6;
      } else {
        SETERRQ(99,"how did I get here?");
      }

    }
  }
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  // communicate b because it will be horizontally differentiated
  ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"IceMISMIPModel: bed topography stored\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceMISMIPModel::additionalAtEndTimestep() {
  // check steady state criterion here
  // write out to data files here  (or just use summaryPrintLine() to do it?)

// FIXME
  return 0;
}


PetscErrorCode IceMISMIPModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

// FIXME:  modify to report grounding line position, av land vel, av shelf vel, g.l. vel.
//         and possibly everything needed for reporting
  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "P         YEAR:     ivol   iarea    meltf     thick0     temp0\n");
    ierr = verbPrintf(2,grid.com,
      "U        years 10^6_km^3 10^6_km^2 (none)          m         K\n");
  } else {
    if (tempAndAge == PETSC_FALSE) {
      ierr = verbPrintf(2,grid.com, "S %12.5f: %8.5f %7.4f   <same> %10.3f    <same>\n",
                         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, H0); CHKERRQ(ierr);
    } else { // general case
      ierr = verbPrintf(2,grid.com, "S %12.5f: %8.5f %7.4f %8.4f %10.3f %9.4f\n",
                         year, volume_kmcube/1.0e6, area_kmsquare/1.0e6, meltfrac,
                         H0,T0); CHKERRQ(ierr);
    }
  }
  return 0;
}

