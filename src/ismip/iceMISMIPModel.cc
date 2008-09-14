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
see run script   examples/mismip/mismip.sh
*/

/* 
WE'VE GOT A PROBLEM WITH FLUX ACROSS GROUNDING LINE; not surprising ...
See "cflx" in output file.

-adapt_ratio is not the whole issue, although with -adapt_ratio 0.12, the default, 
the instability is very noticable as 
highest freq wiggles in dHdt on the land side of the grounding line.  For a Mx=500
case these wiggles extended about 150km back from the grounding line.  Reducing
to "-adapt_ratio 0.08" makes these highest frequency wiggles 
go away, and makes the problem a lower frequency "wobble" very close to the grounding
line (within 40km of it; about 6 to 10 grid spaces).

I think the underlying issue is *regularity of the basal shear stress*, that is, of the 
shear stress coefficient, which jumps from C_MISMIP to zero across the grounding line.
Note that with  export MISMIP_PLAY=1  I can turn on code that makes the ice grounded
all the way out to the calving front, but has zero basal resistance; the results are
similar the intended MISMIP case, even though the surface is not at all what is given
by the floatation criterion.

Damping out C_MISMIP in the 50km on the grounded side of the grounding line (a linear
decrease linearly from C_MISMIP at xg - 50km to zero at xg) does not make a difference,
really.  In other words, I thought we needed
   \tau_b \in W^{1,\infty} 
or something; without the damping all we have is
   \tau_b \in L^\infty.
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


PetscErrorCode MISMIPIce::printInfo(const int thresh, MPI_Comm com) {
  PetscErrorCode ierr;
  ierr = verbPrintf(thresh, com, 
    "Using MISMIP ice w  rho=%6.2f, grav=%6.4f, n=%6.4f, and A=%6.4e.\n",
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


PetscScalar MISMIPIce::softnessParameter(PetscScalar T) const {
  return A_MISMIP;
}


PetscScalar MISMIPIce::hardnessParameter(PetscScalar T) const {
  return B_MISMIP;
}



IceMISMIPModel::IceMISMIPModel(IceGrid &g, IceType *i, MISMIPIce *mismip_i) : 
      IceModel(g, i), mismip_ice(mismip_i) {

  // this just fills the data members with non-junk; see setFromOptions() for decent values
  exper = 1;
  sliding = 'a';
  gridmode = 1;
  runindex = 1;
  runtimeyears = 3.0e4;
  strcpy(initials,"ZZZ");
  m_MISMIP = 1.0/3.0;
  C_MISMIP = 7.624e6;
  regularize_MISMIP = 0.01 / secpera;  
  rstats.xg = -1.0;
}


IceMISMIPModel::~IceMISMIPModel() {
  // close open _t file
  PetscViewerDestroy(tviewfile);
}


PetscErrorCode IceMISMIPModel::printBasalAndIceInfo() {
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
  ierr = mismip_ice->printInfo(2, grid.com); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"IceModel.ice-> info: n=%7.3f, A(T=273.15)=%10.3e\n",
                    ice->exponent(),ice->softnessParameter(273.15)); CHKERRQ(ierr);
  return 0;
}


PetscScalar IceMISMIPModel::basalDragx(PetscScalar **tauc,
                                       PetscScalar **u, PetscScalar **v,
                                       PetscInt i, PetscInt j) const {
  // MAKE SURE THIS IS REALLY BEING USED!!:
  // verbPrintf(1,grid.com,"IceMISMIPModel::basalDragx()\n");
  return basalIsotropicDrag(u, v, i, j);
}


PetscScalar IceMISMIPModel::basalDragy(PetscScalar **tauc,
                                       PetscScalar **u, PetscScalar **v,
                                       PetscInt i, PetscInt j) const {
  return basalIsotropicDrag(u, v, i, j);
}


PetscScalar IceMISMIPModel::basalIsotropicDrag(
            PetscScalar **u, PetscScalar **v, PetscInt i, PetscInt j) const {

  PetscScalar       myC = C_MISMIP;
  
#if (MISMIP_PLAY)
  //verbPrintf(1,grid.com,"doing smoothed C transition");
  const PetscScalar damp_width = 100.0e3;
  if (rstats.xg > damp_width) {  // avoid weird cases where grounding line has moved far left
    // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
    const PetscScalar jfrom0 =
             static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
    const PetscScalar absy = PetscAbs(grid.dy * jfrom0);
    if ((absy > rstats.xg - damp_width) && (absy < rstats.xg)) {
      // this case damps C_MISMIP to zero within damp_width of grounding line
      myC = myC * (rstats.xg - absy) / damp_width;
    }
  }

//  myC = C_MISMIP;  // UNDOES THE ABOVE
#endif

#if (0)
  const PetscScalar switch_to_slick = 710.0e3;
  const PetscScalar jfrom0 =
             static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
  const PetscScalar absy = PetscAbs(grid.dy * jfrom0);
  if (absy > switch_to_slick) {
      //myC = myC / 10.0;
      myC = 0.0;
  }
#endif

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

  // read option    -mismip [1a|1b|2a|2b|3a|3b]
  char Ee[PETSC_MAX_PATH_LEN];
  strcpy(Ee,"");
  ierr = PetscOptionsGetString(PETSC_NULL, "-mismip", Ee, PETSC_MAX_PATH_LEN, PETSC_NULL); 
            CHKERRQ(ierr);
  if (strlen(Ee) != 2) {
    SETERRQ(1,"IceMISMIPModel ERROR:  '-mismip' must be followed by two char argument;\n"
              "  i.e. '-mismip Xx' where Xx=1a,1b,2a,2b,3a,3b");
  } else {
    if ((Ee[0] < '1') || (Ee[0] > '3')) {
      SETERRQ(2,"IceMISMIPModel ERROR:  first character of string 'Xx' in"
                " '-mismip Xx' must be 1, 2, or 3");
    }
    exper = (int) Ee[0] - (int) '0';
    if ((Ee[1] == 'a') || (Ee[1] == 'b')) {
      sliding = Ee[1];
    } else {
      SETERRQ(3,"IceMISMIPModel ERROR:  second character of string 'Xx' in"
                " '-mismip Xx' must be a or b");
    }
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
                                 " should usually be three chars long.");
       CHKERRQ(ierr);
  }

  doTemp                    = PETSC_FALSE;
  doPlasticTill             = PETSC_FALSE;
  doBedDef                  = PETSC_FALSE;

  isDrySimulation           = PETSC_FALSE;
  includeBMRinContinuity    = PETSC_FALSE;
  
  doOceanKill               = PETSC_TRUE;
  
  useSSAVelocity            = PETSC_TRUE;
  computeSurfGradInwardSSA  = PETSC_FALSE;
  useConstantHardnessForSSA = PETSC_FALSE;

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


  // among other things, this reads option "-My"
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  

  // determine gridmode from My
  if (grid.My == 151) 
    gridmode = 1;
  else if (grid.My == 1501) 
    gridmode = 2;
  else
    gridmode = 3;

  // determine model number from doSuperpose; if "-super" then
  //   "ABC2_.." but if not then "ABC1_..."
  if (doSuperpose == PETSC_TRUE) {
    modelnum = 2;
  } else {
    modelnum = 1;
  }
  
  // create prefix (e.g.) "EBU1_2b_M1_A3" for output files with names (e.g.)
  //   EBU1_2b_M1_A3.nc, EBU1_2b_M1_A3_t, EBU1_2b_M1_A3_ss, and EBU1_2b_M1_A3_f
  //   unless user says "-o foo"
  PetscTruth  oused;
  char        oname[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", oname, PETSC_MAX_PATH_LEN, &oused);
           CHKERRQ(ierr);
  if (oused == PETSC_TRUE) {
    strcpy(mprefix, oname);
  } else {
    snprintf(mprefix, sizeof(mprefix), 
             "%s%d_%s_M%d_A%d",
             initials, modelnum, Ee, gridmode, runindex);
    // now act like user set the output name
    ierr = PetscOptionsSetValue("-o",mprefix);  CHKERRQ(ierr);
  }
  ierr = verbPrintf(2,grid.com,
       "IceMISMIPModel:  MISMIP options read.  Will save file %s_t during run, and file\n"
       "  %s.nc at end of run, and files %s_ss, %s_f\n"
       "  at end of run if steady state achieved.\n",
       mprefix,mprefix,mprefix,mprefix); CHKERRQ(ierr);

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
       "IceMISMIPModel: -if or -bif option used; not using"
       "  certain MISMIP formulas to initialize\n");
       CHKERRQ(ierr);
  } else { // usual case: initialize from MISMIP formulas
    ierr = verbPrintf(2,grid.com, 
              "initializing MISMIP experiment %d%c;  grid mode %d;  run %d (A=%5.4e)\n", 
              exper,sliding,gridmode,runindex,mismip_ice->softnessParameter(273.15)); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = createVecs(); CHKERRQ(ierr);

    const PetscScalar   L = 1700.0e3;      // Horizontal half-width of grid
    // NOTE: y takes place of x!!!
    ierr = determineSpacingTypeFromOptions(); CHKERRQ(ierr);

/*    // why does change to Lx=20.0e3 (so cells closer to horizontal aspect ratio of 1)
    //   slow down the linear solve so much?
*/
    // is Lz = 4000m an adequate choice for thickness for all runs?
    //   (could be set in setFromOptions() according to experiment/run/...)
    // effect of double rescale is to compute grid.dy so we can get square cells
    //   (in horizontal)
    ierr = grid.rescale_and_set_zlevels(1000.0e3, L, 4000.0,PETSC_TRUE,PETSC_FALSE); CHKERRQ(ierr); 
    const PetscScalar   Lx_desired = (grid.dy * grid.Mx) / 2.0;
    ierr = grid.rescale_and_set_zlevels(Lx_desired, L, 4000.0,PETSC_TRUE,PETSC_FALSE); CHKERRQ(ierr); 

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

    initialized_p = PETSC_TRUE;
  }
  
  ierr = IceModel::initFromOptions(PETSC_TRUE); CHKERRQ(ierr);  // regridding can happen here

  if (!isInitialized()) {
    SETERRQ(1, "ERROR: IceMISMIPModel has not been initialized!\n");
  }

  // use MISMIP runtimeyears UNLESS USER SPECIFIES A RUN LENGTH
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

  ierr = printBasalAndIceInfo(); CHKERRQ(ierr);
  
// report on these flags: doTemp=false, doBedDef=false, doPlasticTill=false  ?
// check "-ssa" is set? check on -super option?

  // create ABC1_1b_M1_A1_t kind of file for every 50 year results
  strcpy(tfilename,mprefix);
  strcat(tfilename,"_t");
  ierr = PetscViewerASCIIOpen(grid.com, tfilename, &tviewfile); CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(tviewfile, PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);

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
#if (0)
        if (PetscAbs(y) < 710.0e3) {
          b[i][j] = 720.0 - 778.5 * xs;
        } else {
          b[i][j] = 720.0 - 778.5 * (710.0e3 / 750.0e3);
        } 
#else
        b[i][j] = 720.0 - 778.5 * xs;
#endif
      } else if (exper == 3) {
        const PetscScalar xs2 = xs * xs,
                          xs4 = xs2 * xs2,
                          xs6 = xs4 * xs2;
        b[i][j] = 729.0 - 2184.0 * xs2 + 1031.72 * xs4 - 151.72 * xs6;
      } else {
        SETERRQ(99,"how did I get here?");
      }

#if (0)
      b[i][j] += 940.8; // eliminate the shelf by putting the calving front right where bed=0
                        // (note mask = MASK_FLOATING_OCEAN0 beyond calving front)
                        // calculation was   A + 720 - 778.5 * (1600e3/750e3) = 0 
                        // to get A = 940.8
#endif
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

  const PetscScalar MISMIP_calving_front = 1600.0e3;

  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
    
      // NOTE !!!!:   y  REPLACES   x   FOR VIEWING CONVENIENCE!
      const PetscScalar jfrom0 =
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
      const PetscScalar y = grid.dy * jfrom0;
      if (PetscAbs(y) >= MISMIP_calving_front) {
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

  PetscErrorCode  ierr;

  // why is it necessary to globalize this norm??  because vdHdt is local??
  PetscScalar     infnormdHdt, ginfnormdHdt;
  ierr = VecNorm(vdHdt,NORM_INFINITY,&infnormdHdt); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&infnormdHdt, &ginfnormdHdt, grid.com); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"max|dH/dt| = %8.3e m/a; MISMIP steady is < 1.0e-4 m/a");
     CHKERRQ(ierr);

  if (ginfnormdHdt < 1.0e-4 / secpera) {  // if all points have dHdt < 10^-4 m/yr,
    // then set the IceModel goal of endYear to the current year; stops immediately
    endYear = grid.year;  

    // report stopping to standard out
    ierr = verbPrintf(2,grid.com,
        "\nIceMISMIPModel: MISMIP steady state criterion (max|dH/dt| < 10^-4 m/yr) satisfied;\n"
        "                stopping at year=%.3f\n",grid.year); CHKERRQ(ierr);

    // leave stopping stamp in output NetCDF file
    char str[HISTORY_STRING_LENGTH];
    snprintf(str, sizeof(str), 
       "MISMIP steady state criterion (max|dHdt| < 10^-4 m/yr) satisfied.  Stopping."
       "  Completed timestep at year=%.3f.",grid.year);
    stampHistory(str); 

    // get stats in preparation for writing
    ierr = getRoutineStats(); CHKERRQ(ierr);
    ierr = getMISMIPStats(); CHKERRQ(ierr);

    // write ASCII file ABC1_1b_M1_A1_ss
    PetscViewer  viewfile;
    char         filename[PETSC_MAX_PATH_LEN];
    strcpy(filename,mprefix);
    strcat(filename,"_ss");
    ierr = verbPrintf(2, grid.com, "IceMISMIPModel:  writing profile to %s\n",
             filename); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(grid.com, filename, &viewfile); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewfile, PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);

    PetscScalar     **H;
    ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
    if (grid.xs == 0) {  // if (0,jg) is in ownership
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        const PetscScalar jfrom0 =
               static_cast<PetscScalar>(j)-static_cast<PetscScalar>(grid.My - 1)/2.0;
        const PetscScalar y = grid.dy * jfrom0;
        if (y >= 0) {
          ierr = PetscViewerASCIISynchronizedPrintf(viewfile,
                 "%10.2f %10.4f\n", y / 1000.0, H[0][j]); CHKERRQ(ierr);
        }
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
    
    ierr = PetscViewerFlush(viewfile); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewfile); CHKERRQ(ierr);
    
    // write ASCII file ABC1_1b_M1_A1_f
    strcpy(filename,mprefix);
    strcat(filename,"_f");
    ierr = verbPrintf(2, grid.com, 
             "IceMISMIPModel:  writing final location of grounding line to %s\n",
             filename); CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(grid.com, filename, &viewfile); CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewfile, PETSC_VIEWER_ASCII_DEFAULT); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewfile,"%10.4f %10.2f\n", rstats.xg / 1000.0, grid.year);
               CHKERRQ(ierr);
    ierr = PetscViewerDestroy(viewfile); CHKERRQ(ierr);

  }
  return 0;
}


PetscErrorCode IceMISMIPModel::getMISMIPStats() {
  // run this only after getRoutineStats() is called
  
  PetscErrorCode  ierr;
  PetscScalar     **H, **b, **q;
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  ierr = VecPointwiseMult(vWork2d[0], vvbar, vH); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &q); CHKERRQ(ierr);
  // q[i][j] is signed flux in y direction, in units of m^2/s
  
  mstats.x1 = rstats.xg;
  mstats.x2 = rstats.xg - grid.dy;
  mstats.x3 = rstats.xg + grid.dy;
  
  PetscScalar myh2 = 0.0, myh3 = 0.0, myb1 = -1e6, myb2 = -1e6, myb3 = -1e6, 
              myq1 = -1e20, myq2 = -1e20, myq3 = -1e20;
  const int jg = (int)floor(rstats.jg + 0.1);

  mstats.h1 = rstats.hxg;  // already computed
  if ( (jg >= grid.ys) && (jg < grid.ys + grid.ym)
       && (grid.xs == 0)                             ) {  // if (0,jg) is in ownership
    myb1 = b[0][jg];
    myq1 = q[0][jg];
  }
  ierr = PetscGlobalMax(&myb1, &mstats.b1, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&myq1, &mstats.q1, grid.com); CHKERRQ(ierr);

  if ( (jg-1 >= grid.ys) && (jg-1 < grid.ys + grid.ym)
       && (grid.xs == 0)                             ) {  // if (0,jg-1) is in ownership
    myh2 = H[0][jg-1];
    myb2 = b[0][jg-1];
    myq2 = q[0][jg-1];
  }
  ierr = PetscGlobalMax(&myh2, &mstats.h2, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&myb2, &mstats.b2, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&myq2, &mstats.q2, grid.com); CHKERRQ(ierr);

  if ( (jg+1 >= grid.ys) && (jg+1 < grid.ys + grid.ym)
       && (grid.xs == 0)                             ) {  // if (0,jg+1) is in ownership
    myh3 = H[0][jg+1];
    myb3 = b[0][jg+1];
    myq3 = q[0][jg+1];
  }
  ierr = PetscGlobalMax(&myh3, &mstats.h3, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&myb3, &mstats.b3, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalMax(&myq3, &mstats.q3, grid.com); CHKERRQ(ierr);

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &q); CHKERRQ(ierr);

  // perform MISMIP diagnostic computation here, to estimate dxg/dt:
  //   d xg            a - dq/dx
  //   ---- = -----------------------------
  //    dt     dh/dx - (rhow/rhoi) (db/dx)
  const PetscScalar dqdx = (mstats.q1 - mstats.q2) / (mstats.x1 - mstats.x2),
                    dhdx = (mstats.h1 - mstats.h2) / (mstats.x1 - mstats.x2),
                    dbdx = (mstats.b1 - mstats.b2) / (mstats.x1 - mstats.x2);
  mstats.dxgdt = ((0.3/secpera) - dqdx) / (dhdx - (ocean.rho/mismip_ice->rho) * dbdx);  
  return 0;
}


PetscErrorCode IceMISMIPModel::getRoutineStats() {
  PetscErrorCode  ierr;

  PetscScalar     **mask, **H, **vbar;

  // these are in MKS; sans "g" are local to the processor; with "g" are global 
  //   across all processors; we only evaluate for x > 0
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
  rstats.jg = gjg;
  
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

/*
Because this model does resolve the shelf and only uses the floatation criterion
to move the grounding line, we will give 17 numbers.  These numbers will go into
a reportable ascii file ABC1_1a_M1_A1_t.

The reported numbers are

   t  x_g  Volume  h(0,t)  h(x_g,t)
     x_1 h(x_1,t) b(x_1) q(x_1,t)       // last grounded point (i.e. x_1 = x_g)
     x_2 h(x_2,t) b(x_2) q(x_2,t)       // x_2 = x_1 - dx
     x_3 h(x_3,t) b(x_3) q(x_3,t)       // x_3 = x_1 + dx

The number of reported digits is

     8 chars (includes ".") for t [years]
     7 chars for x_*       [km]
     8 chars for Volume    [10^6 km^3]
     7 chars for h(*,t)    [m]
     7 chars for b(*,t)    [m]
     7 chars for q(*,t)    [m^2/year]

The format written to ABC1_1a_M1_A1_t has 137 columns, like this:

######## ####### ######## ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

The at verbosity 3 or higher, the stdout format fits in an 80 column line:

M  ######## ####### ######## ####### #######
   ####### ####### ####### ####### ####### ####### ####### #######
   ####### ####### ####### #######

A line
   [ d(xg)/dt = ####### m/yr by MISMIP computation ]
is written to stdout at verbosity 2 or higher.  This grounding line motion rate
is computed as in MISMIP description, and finite differences.
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
    if (fabs(fmod(year, 50.0)) < 1.0e-6) {
      // write another line to ASCII file ABC1_1b_M1_A1_t; also write to stdout
      //   (given some verbosity level); note modest code redundancy
      ierr = verbPrintf(2, grid.com, 
             "[IceMISMIPModel:  adding t=%10.3f line to file %s;\n",
             year,tfilename); CHKERRQ(ierr);
      ierr = getMISMIPStats(); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com,"M  ");
      ierr = verbPrintf(3,grid.com,
        "%8.2f %7.2f %8.5f %7.2f %7.2f ",
        year, rstats.xg / 1000.0, volume_kmcube/1.0e6, H0, rstats.hxg); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tviewfile,
        "%8.2f %7.2f %8.5f %7.2f %7.2f ",
        year, rstats.xg / 1000.0, volume_kmcube/1.0e6, H0, rstats.hxg); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com,"\n   ");
      ierr = verbPrintf(3,grid.com,
        "%7.2f %7.2f %7.2f %7.0f %7.2f %7.2f %7.2f %7.0f ",
        mstats.x1 / 1000.0, mstats.h1, mstats.b1, mstats.q1 * secpera,
        mstats.x2 / 1000.0, mstats.h2, mstats.b2, mstats.q2 * secpera); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tviewfile,
        "%7.2f %7.2f %7.2f %7.0f %7.2f %7.2f %7.2f %7.0f ",
        mstats.x1 / 1000.0, mstats.h1, mstats.b1, mstats.q1 * secpera,
        mstats.x2 / 1000.0, mstats.h2, mstats.b2, mstats.q2 * secpera); CHKERRQ(ierr);
      ierr = verbPrintf(3,grid.com,"\n   ");
      ierr = verbPrintf(3,grid.com,
        "%7.2f %7.2f %7.2f %7.0f\n",
        mstats.x3 / 1000.0, mstats.h3, mstats.b3, mstats.q3 * secpera); CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(tviewfile,
        "%7.2f %7.2f %7.2f %7.0f\n",
        mstats.x3 / 1000.0, mstats.h3, mstats.b3, mstats.q3 * secpera); CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com,
        "   d(xg)/dt = %10.2f m/yr by MISMIP computation ]\n",
        mstats.dxgdt * secpera); CHKERRQ(ierr);
      
    }
  }
  return 0;
}

