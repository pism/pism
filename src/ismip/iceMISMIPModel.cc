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
WORKING REFERENCE RUNs:

  mpiexec -n 2 pisms -mismip 1b -Mx 5 -My 101 -Mz 11 -super -ocean_kill \
     -ksp_rtol 1e-11 -ssa_rtol 1e-6 -y 10000

and version without -super, and with higher res
*/

/* 

WE'VE GOT A PROBLEM WITH FLUX ACROSS GROUNDING LINE; not surprising ...
See "cflx" in output file.

-adapt_ratio is not the whole issue, although the instability is very noticable as 
highest freq wiggles in dHdt on the land side of the grounding line.  For a Mx=500
case these wiggles extended about 150km back from the grounding line.  Reducing
to "-adapt_ratio 0.08" (from 0.12, default) makes these highest frequency wiggles 
go away, and makes the problem a lower frequency "wobble" very close to the grounding
line (within 40km of it; about 6 to 10 grid spaces).

I think the underlying issue is *regularity of the basal shear stress*, that is, of the 
shear stress coefficient, which jumps from C_MISMIP to zero across the grounding line.

A plan: in additionalAtEndTimestep() we set an IceMISMIPModel variable which is the 
current estimate of the location of the grounding line.  It must be determined correctly
in parallel!  In particular, it will be xg = x_{i+1/2} for the largest i for which mask(i)
is DRAGGING *and* mask(i+1) is FLOATING.  If this impossible it will be negative.
PetscGlobalMax() should be usable.  Then we will damp out C_MISMIP in the 50km on the 
grounded side of the grounding line; it will decrease linearly from C_MISMIP at xg - 50km
to zero at xg.

In other words, I think we need
   \tau_b \in W^{1,\infty} 
or something; without the damping all we have is
   \tau_b \in L^\infty.
   
*/

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
    "Using MISMIP ice w  rho=%5.2f, grav=%5.4f, n=%d, and A=%5.4e.\n",
    rho, grav, n, A_MISMIP); CHKERRQ(ierr);
  return 0;
}


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


IceMISMIPModel::IceMISMIPModel(IceGrid &g, IceType *i, MISMIPIce *mismip_i) : 
      IceModel(g, i), mismip_ice(mismip_i) {

  // this just fills the data members with non-junk; see setFromOptions() for decent values
  exper = 1;
  sliding = 'a';
  gridmode = 1;
  runindex = 1;
  runtimeyears = 3.0e4;
  strcpy(initials,"ABC");
  m_MISMIP = 1.0/3.0;
  C_MISMIP = 7.624e6;
  regularize_MISMIP = 0.01 / secpera;  
  rstats.xg = -1.0;
}


PetscErrorCode IceMISMIPModel::printBasalInfo() {
  PetscErrorCode ierr;
  if (m_MISMIP == 1.0) {
    ierr = verbPrintf(2, grid.com, 
      "Using MISMIP sliding w  tau_b = - C u,  C=%5.4e.\n", C_MISMIP); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, 
      "Using MISMIP sliding w  tau_b = - C (|u|^2 + eps^2)^{(m-1)/2} u,\n"
      "   m=%5.4f, C=%5.4e, and eps = %5.4f m/a.\n",
      m_MISMIP, C_MISMIP, regularize_MISMIP * secpera); CHKERRQ(ierr);
  }
  //ierr = verbPrintf(2, grid.com, 
  // "[printBasalInfo: rstats.xg=%5.4f]\n", rstats.xg); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceMISMIPModel::basalDragx(PetscScalar **beta, PetscScalar **tauc,
                                       PetscScalar **u, PetscScalar **v,
                                       PetscInt i, PetscInt j) const {
  return basalIsotropicDrag(beta, tauc, u, v, i, j);
}


PetscScalar IceMISMIPModel::basalDragy(PetscScalar **beta, PetscScalar **tauc,
                                       PetscScalar **u, PetscScalar **v,
                                       PetscInt i, PetscInt j) const {
  return basalIsotropicDrag(beta, tauc, u, v, i, j);
}


PetscScalar IceMISMIPModel::basalIsotropicDrag(PetscScalar **beta, PetscScalar **tauc,
                                       PetscScalar **u, PetscScalar **v,
                                       PetscInt i, PetscInt j) const {
  PetscScalar myC = C_MISMIP;
  if (rstats.xg < 50.0) {
    // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
    const PetscScalar jfrom0 =
             static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
    const PetscScalar absy = PetscAbs(grid.dy * jfrom0);
    if ((absy > rstats.xg - 50.0e3) && (absy < rstats.xg)) {
      myC = myC * (rstats.xg - absy) / 50.0e3;
    }
  }

  if (m_MISMIP == 1.0) {
    return myC;
  } else {
    const PetscScalar magsliding = PetscSqr(u[i][j]) + PetscSqr(v[i][j])
                                   + PetscSqr(regularize_MISMIP);
    return myC * pow(magsliding, (m_MISMIP - 1.0) / 2.0);
  }
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
      SETERRQ(1,"IceMISMIPModel ERROR:  first character of string Ee in"
                " '-mismip Ee' must be 1, 2, or 3");
    }
    exper = (int) Ee[0] - (int) '0';
    sliding = 'a';
  } else if (strlen(Ee) == 2) {
    if ((Ee[0] < '1') || (Ee[0] > '3')) {
      SETERRQ(1,"IceMISMIPModel ERROR:  first character of string Ee in"
                " '-mismip Ee' must be 1, 2, or 3");
    }
    exper = (int) Ee[0] - (int) '0';
    if ((Ee[1] == 'a') || (Ee[1] == 'b')) {
      sliding = Ee[1];
    } else {
      SETERRQ(2,"IceMISMIPModel ERROR:  second character of string Ee in"
                " '-mismip Ee' must be a or b");
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
      SETERRQ(5,"IceMISMIPModel ERROR:  run index N in '-run N' must be"
                " <= 9 in experiments 1 or 2");
    }
    runtimeyears = 3.0e4;
    ierr = mismip_ice->setA(Aexper1or2[runindex]); CHKERRQ(ierr);  
  } else if (exper == 3) {
    if (sliding == 'a') {
      if (runindex > 13) {
        SETERRQ(6,"IceMISMIPModel ERROR:  run index N in '-run N' must be"
                  " <= 13 in experiment 3a");
      }
      runtimeyears = timeexper3a[runindex];
      ierr = mismip_ice->setA(Aexper3a[runindex]); CHKERRQ(ierr);  
    } else if (sliding == 'b') {
      if (runindex > 15) {
        SETERRQ(7,"IceMISMIPModel ERROR:  run index N in '-run N' must be"
                  " <= 15 in experiment 3b");
      }
      runtimeyears = timeexper3b[runindex];
      ierr = mismip_ice->setA(Aexper3b[runindex]); CHKERRQ(ierr);  
    } else {
      SETERRQ(99, "how did I get here?");
    }
  }

  // read option  -initials ABC
  ierr = PetscOptionsGetString(PETSC_NULL, "-initials", initials, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (strlen(initials) != 3) {
    ierr = verbPrintf(1,grid.com,"IceMISMIPModel WARNING:  Initials string"
                                 " should be three chars long.");
       CHKERRQ(ierr);
  }

  doTemp               = PETSC_FALSE;
  doPlasticTill        = PETSC_FALSE;
  doBedDef             = PETSC_FALSE;

  isDrySimulation        = PETSC_FALSE;
  includeBMRinContinuity = PETSC_FALSE;
  
  useSSAVelocity            = PETSC_TRUE;
  computeSurfGradInwardSSA  = PETSC_TRUE;
  useConstantHardnessForSSA = PETSC_FALSE;

// use value of gridmode to set grid??  set grid.Mx, grid.My from options or lack thereof

// check "-ssa" is set

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
       "IceMISMIPModel WARNING: -if or -bif option used; not using"
       " MISMIP formulae to initialize\n");
       CHKERRQ(ierr);  
  } else { // usual case: initialize from MISMIP formulas
    ierr = verbPrintf(1,grid.com, 
              "initializing MISMIP experiment %d%c\n"
              "    [grid mode %d; run %d (A=%5.4e, MISMIP run length = %5.2f yrs)]\n", 
              exper,sliding,gridmode,runindex,mismip_ice->getA(),runtimeyears); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);

    // see Table 3
    if (sliding == 'a') {
      m_MISMIP = 1.0/3.0;
      C_MISMIP = 7.624e6;
    } else if (sliding == 'b') {
      m_MISMIP = 1.0;
      C_MISMIP = 7.2082e10;
    } else {
      SETERRQ(99, "how did I get here?");
    }
    regularize_MISMIP = 0.01 / secpera;  

    // FIXME:  use value of gridmode to set grid??

    const PetscScalar   L = 1700.0e3;      // Horizontal half-width of grid
    // NOTE: y takes place of x!!!
    // NOTE: put calving front at 1600.0e3  ??
    ierr = determineSpacingTypeFromOptions(); CHKERRQ(ierr);

    // FIXME:  what to put for thickness here?
    ierr = grid.rescale_and_set_zlevels(200.0e3, L, 4000.0); CHKERRQ(ierr); 

    // all of these relate to models which should be turned off ...
    ierr = VecSet(vHmelt, 0.0); CHKERRQ(ierr);
    // none use Goldsby-Kohlstedt or do age calc
    setInitialAgeYears(initial_age_years_default);
    ierr = VecSet(vuplift,0.0); CHKERRQ(ierr);  // no bed deformation
    ierr = VecSet(vTs, ice->meltingTemp); CHKERRQ(ierr);
    ierr = T3.setToConstant(ice->meltingTemp); CHKERRQ(ierr);
    ierr = Tb3.setToConstant(ice->meltingTemp); CHKERRQ(ierr);

    ierr = VecSet(vAccum, 0.3/secpera); CHKERRQ(ierr);

    ierr = VecSet(vH, 10.0); CHKERRQ(ierr);  // initial thickness of 10 m

    ierr = setMISMIPBed(); CHKERRQ(ierr);
    ierr = setMISMIPMask(); CHKERRQ(ierr);
    ierr = verbPrintf(4,grid.com,"IceMISMIPModel: bed topography and mask stored\n");
              CHKERRQ(ierr);

    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

    //    WHEN INITIALIZING W/O INPUT FILE, USE MISMIP runtimeyears UNLESS USER
    //      SPECIFIES A RUN LENGTH
    // use -y option, if given, to overwrite runtimeyears
    PetscTruth ySet, ysSet, yeSet;
    ierr = PetscOptionsHasName(PETSC_NULL, "-y", &ySet); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-ys", &ysSet); CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL, "-ye", &yeSet); CHKERRQ(ierr);
    if ( (ySet == PETSC_TRUE) || ( (ysSet == PETSC_TRUE) && (yeSet == PETSC_TRUE) ) ) {
      ierr = verbPrintf(2,grid.com,
        "IceMISMIPModel: ignoring MISMIP-specified run length and using value\n"
        "  from user option -y (or -ys and -ye)\n"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com,
        "IceMISMIPModel: setting run length to %5.2f years (from MISMIP specs)\n",
        runtimeyears); CHKERRQ(ierr);
      grid.year = 0.0;
      ierr = setStartYear(grid.year); CHKERRQ(ierr);
      ierr = setEndYear(grid.year + runtimeyears); CHKERRQ(ierr);
      yearsStartRunEndDetermined = PETSC_TRUE;
    }

    initialized_p = PETSC_TRUE;
  }
  
  ierr = mismip_ice->printInfo(2,grid.com); CHKERRQ(ierr);
//  ierr = basal->printInfo(2,grid.com); CHKERRQ(ierr);
  ierr = printBasalInfo(); CHKERRQ(ierr);
  
  if (!isInitialized()) {
    SETERRQ(1, "ERROR: IceMISMIPModel has not been initialized!\n");
  }

  ierr = IceModel::initFromOptions(PETSC_TRUE); CHKERRQ(ierr);  // regridding can happen here

  // FIXME: report on these flags: doTemp=false, doBedDef=false, doPlasticTill=false

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
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
      const PetscScalar y = grid.dy * jfrom0;
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
  return 0;
}



PetscErrorCode IceMISMIPModel::setMISMIPMask() {
  PetscErrorCode ierr;
  PetscScalar    **mask;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
      const PetscScalar jfrom0 =
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
      const PetscScalar y = grid.dy * jfrom0;
      if (PetscAbs(y) >= 1600.0e3) {
        mask[i][j] = MASK_FLOATING_OCEAN0;
      } else {
        // note updateSurfaceElevationAndMask() will mark DRAGGING as FLOATING if it is floating
        mask[i][j] = MASK_DRAGGING;
      }

    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  // communicate it
  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceMISMIPModel::additionalAtStartTimestep() {
  // this is called at the beginning of each pass through time-stepping loop in IceModel::run()

  // go to next multiple of 50 years
  const PetscScalar tonext50 = (50.0 - fmod(grid.year, 50.0)) * secpera;
  if (maxdt_temporary < 0.0) {  // it has not been set
    maxdt_temporary = tonext50;
  } else {
    maxdt_temporary = PetscMin(maxdt_temporary, tonext50);
  }
  return 0;
}


PetscErrorCode IceMISMIPModel::additionalAtEndTimestep() {
  // this is called at the end of each pass through time-stepping loop in IceModel::run()

// FIXME: check steady state criterion here
  return 0;
}


PetscErrorCode IceMISMIPModel::getMISMIPStats() {
  // PetscErrorCode  ierr;

// FIXME: see summaryPrintLine() for what is needed
  return 0;
}


PetscErrorCode IceMISMIPModel::getRoutineStats() {
  PetscErrorCode  ierr;

  PetscScalar     **mask, **H, **vbar;

  // these are in MKS; sans "g" are local to the processor; with "g" are global across all processors
  // we only evaluate for x > 0
  PetscScalar     maxubar = 0.0, avubargrounded = 0.0, avubarfloating = 0.0, jg = 0.0,
                  Ngrounded = 0.0, Nfloating = 0.0;
  PetscScalar     gavubargrounded, gavubarfloating, gjg,
                  gNgrounded, gNfloating;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      const PetscScalar jfrom0 =
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;

      // grounding line xg is largest  y  so that  mask[i][j] != FLOATING
      //       and mask[i][j+1] == FLOATING
      if ( (jfrom0 > 0.0) && (H[i][j] > 0.0) 
           && (modMask(mask[i][j]) != MASK_FLOATING) 
           && (modMask(mask[i][j+1]) == MASK_FLOATING) ) {
        // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
        jg = PetscMax(jg,static_cast<PetscScalar>(j));
      }

      if ((jfrom0 > 0) && (H[i][j] > 0.0)) {
        // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
        if (vbar[i][j] > maxubar)  maxubar = vbar[i][j];
        if (modMask(mask[i][j]) != MASK_FLOATING) {
          Ngrounded += 1.0;
          avubargrounded += vbar[i][j];
        } else {
          Nfloating += 1.0;
          avubarfloating += vbar[i][j];        
        }
      }

    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbar, &vbar); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&jg, &gjg, grid.com); CHKERRQ(ierr);
  const PetscScalar gjgfrom0 =
          gjg - static_cast<PetscScalar>(grid.My - 1)/2.0;
  rstats.xg = gjgfrom0 * grid.dy;
  
  PetscScalar myhxg = 0.0;
  if ( (gjg >= grid.ys) && (gjg < grid.ys + grid.ym)
       && (grid.xs == 0)                             ) {  // if (0,gjg) is in ownership
    myhxg = H[0][static_cast<int>(gjg)]; // i.e. hxg = H[0][gjg]
  } else {
    myhxg = 0.0;
  } 
  ierr = PetscGlobalMax(&myhxg, &rstats.hxg, grid.com); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&maxubar, &rstats.maxubar, grid.com); CHKERRQ(ierr);
    
  ierr = PetscGlobalSum(&Ngrounded, &gNgrounded, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avubargrounded, &gavubargrounded, grid.com); CHKERRQ(ierr);
  if (gNgrounded > 0)   gavubargrounded = gavubargrounded / gNgrounded;
  else                  gavubargrounded = 0.0;  // degenerate case
  rstats.avubarG = gavubargrounded;

  ierr = PetscGlobalSum(&Nfloating, &gNfloating, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&avubarfloating, &gavubarfloating, grid.com); CHKERRQ(ierr);
  if (gNfloating > 0)   gavubarfloating = gavubarfloating / gNfloating;
  else                  gavubarfloating = 0.0;  // degenerate case
  rstats.avubarF = gavubarfloating;

  return 0;
}

PetscErrorCode IceMISMIPModel::summaryPrintLine(
    const PetscTruth printPrototype, const PetscTruth tempAndAge,
    const PetscScalar year, const PetscScalar dt, 
    const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
    const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0) {

// FIXME:  modify to report everything needed for reporting, with leading "M "

/*
Because this model does resolve the shelf and only uses the floatation criterion
to move the grounding line, we will give 17 numbers.
(These numbers will go into ABC1_1a_M1_A1_t.)

The basic format is

   t  x_g  Volume  h(0,t)  h(x_g,t)
     x_1 h(x_1,t) b(x_1) q(x_1,t)       // last grounded point (i.e. x_1 = x_g)
     x_2 h(x_2,t) b(x_2) q(x_2,t)       // x_2 = x_1 - dx
     x_3 h(x_3,t) b(x_3) q(x_3,t)       // x_3 = x_1 + dx

The actual format is

M  ######## ####### ######## ####### #######
    ####### ####### ####### #######
    ####### ####### ####### #######
    ####### ####### ####### #######

with 8 chars (includes ".") for t [years]
     7 chars for x_*       [km]
     8 chars for Volume    [10^6 km^3]
     7 chars for h(*,t)    [m]
     7 chars for b(*,t)    [m]
     7 chars for q(*,t)    [m^2/year]
*/

  PetscErrorCode ierr;
  if (printPrototype == PETSC_TRUE) {
    ierr = verbPrintf(2,grid.com,
      "P         YEAR:     ivol      h0      xg     hxg maxubar avubarG avubarF\n");
    ierr = verbPrintf(2,grid.com,
      "U        years 10^6_km^3       m      km       m     m/a     m/a     m/a\n");
  } else {
    ierr = getRoutineStats(); CHKERRQ(ierr);
    ierr = verbPrintf(2,grid.com,
      "S %12.5f: %8.5f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
      year, volume_kmcube/1.0e6, 
      H0, rstats.xg / 1000.0, rstats.hxg, rstats.maxubar * secpera, 
      rstats.avubarG * secpera, rstats.avubarF * secpera); CHKERRQ(ierr);
    /*  
    if (CRITERIA FOR HITTING mult. of 50 year)
      ierr = getMISMIPStats(); CHKERRQ(ierr);
      PRINT "M " summary
    */
  }

  return 0;
}

