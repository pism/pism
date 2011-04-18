// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "grid.hh"
#include "iceModel.hh"
#include "iceEISModel.hh"
#include "SIAFD.hh"
#include "SIA_Sliding.hh"

IceEISModel::IceEISModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides)
  : IceModel(g, conf, conf_overrides) {
  expername = 'A';

  // following flag must be here in constructor because IceModel::createVecs()
  //   uses it; can't wait till init_physics()
  // non-polythermal methods; can be overridden by the command-line option -no_cold:
  config.set_flag("do_cold_ice_methods", true);
}

PetscErrorCode IceEISModel::createVecs() {
  PetscErrorCode ierr = IceModel::createVecs(); CHKERRQ(ierr);

  // this ensures that these variables are saved to an output file and are read
  // back in if -i option is used (they are "model_state", in a sense, since
  // PSDummy is used):
  ierr = artm.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  ierr = acab.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  ierr = variables.add(acab); CHKERRQ(ierr); 
  ierr = variables.add(artm); CHKERRQ(ierr); 

  return 0;
}

//! Only executed if NOT initialized from file (-i).
PetscErrorCode IceEISModel::set_grid_defaults() {
  grid.Lx = 750e3;
  grid.Ly = 750e3;
  grid.Lz = 4e3;  // depend on auto-expansion to handle bigger thickness
  return 0;
}


//! Option -eisII determines the single character name of EISMINT II experiments.
/*! Example is "-eisII F".   Defaults to experiment A.  */
PetscErrorCode IceEISModel::set_expername_from_options() {
  PetscErrorCode      ierr;

  string eisIIexpername = "A";
  char temp = expername;
  bool EISIIchosen;
  ierr = PISMOptionsString("-eisII", "EISMINT II experiment name",
			   eisIIexpername, EISIIchosen);
            CHKERRQ(ierr);
  if (EISIIchosen == PETSC_TRUE) {
    temp = (char)toupper(eisIIexpername.c_str()[0]);
    if ((temp >= 'A') && (temp <= 'L')) {
      expername = temp;
    } else {
      ierr = PetscPrintf(grid.com,
        "option -eisII must have value A, B, C, D, E, F, G, H, I, J, K, or L\n");
        CHKERRQ(ierr);
      PISMEnd();
    }
  }

  char tempstr[TEMPORARY_STRING_LENGTH];
  snprintf(tempstr, TEMPORARY_STRING_LENGTH, "%c", temp);

  config.set_string("EISMINT_II_experiment", tempstr);

  return 0;
}


PetscErrorCode IceEISModel::setFromOptions() {
  PetscErrorCode      ierr;

  ierr = PetscOptionsBegin(grid.com, "", "IceEISModel options", ""); CHKERRQ(ierr);

  ierr = set_expername_from_options(); CHKERRQ(ierr);

  config.set("enhancement_factor", 1.0);
  config.set("bed_smoother_range", 0.0);  // none use bed smoothing & bed roughness
                                          // parameterization
  // basal melt does not change computation of mass continuity or vertical velocity:
  config.set_flag("include_bmr_in_continuity", false);

  ierr = verbPrintf(2,grid.com, 
    "setting parameters for surface mass balance and temperature in EISMINT II experiment %c ... \n", 
    expername); CHKERRQ(ierr);
  // EISMINT II specified values for parameters
  S_b = 1.0e-2 * 1e-3 / secpera;    // Grad of accum rate change
  S_T = 1.67e-2 * 1e-3;           // K/m  Temp gradient
  // these are for A,E,G,H,I,K:
  M_max = 0.5 / secpera;  // Max accumulation
  R_el = 450.0e3;           // Distance to equil line (accum=0)
  T_min = 238.15;
  switch (expername) {
    case 'B':  // supposed to start from end of experiment A and:
      T_min = 243.15;
      break;
    case 'C':
    case 'J':
    case 'L':  // supposed to start from end of experiment A (for C;
               //   resp I and K for J and L) and:
      M_max = 0.25 / secpera;
      R_el = 425.0e3;
      break;
    case 'D':  // supposed to start from end of experiment A and:
      R_el = 425.0e3;
      break;
    case 'F':  // start with zero ice and:
      T_min = 223.15;
      break;
  }

  // if user specifies Tmin, Tmax, Mmax, Sb, ST, Rel, then use that (override above)
  bool paramSet;
  ierr = PISMOptionsReal("-Tmin", "T min, Kelvin",
			 T_min, paramSet); CHKERRQ(ierr);
  ierr = PISMOptionsReal("-Tmax", "T max, Kelvin",
			 T_max, paramSet); CHKERRQ(ierr);

  PetscReal myMmax = M_max*secpera,
    mySb = S_b * secpera / 1e3,
    myST = S_T / 1e3,
    myRel = R_el / 1e3;
  ierr = PISMOptionsReal("-Mmax", "Maximum accumulation, m/year",
			 myMmax, paramSet); CHKERRQ(ierr);
  if (paramSet)     M_max = myMmax / secpera;

  ierr = PISMOptionsReal("-Sb", "FIXME",
			 mySb, paramSet); CHKERRQ(ierr);
  if (paramSet)     S_b = mySb * 1e-3 / secpera;

  ierr = PISMOptionsReal("-ST", "FIXME",
			 myST, paramSet); CHKERRQ(ierr);
  if (paramSet)     S_T = myST * 1e-3;

  ierr = PISMOptionsReal("-Rel", "km; FIXME",
			 myRel, paramSet); CHKERRQ(ierr);
  if (paramSet)     R_el = myRel * 1e3;

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // these change default IceModel semantics to match EISMINT II choices:

  bool flag;
  ierr = PISMOptionsIsSet("-track_Hmelt", flag); CHKERRQ(ierr);
  if (flag) updateHmelt = PETSC_TRUE;

  ierr = IceModel::setFromOptions();  CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEISModel::init_physics() {
  PetscErrorCode ierr;

  iceFactory.setType(ICE_PB);  // Paterson-Budd; can be overridden by options

  if (ice == NULL) {
    // Initialize the IceFlowLaw object:
    if (config.get_flag("do_cold_ice_methods") == false) {
      ierr = verbPrintf(2, grid.com,
                        "  setting flow law to polythermal type ...\n"); CHKERRQ(ierr);
      ierr = verbPrintf(3, grid.com,
                        "      (= Glen-Paterson-Budd-Lliboutry-Duval type)\n"); CHKERRQ(ierr);

      // new flowlaw which has dependence on enthalpy, not temperature
      ice = new GPBLDIce(grid.com, "", config);

    } else {
      ierr = verbPrintf(2, grid.com,
                        "  doing cold ice methods ...\n"); CHKERRQ(ierr);

      ierr = iceFactory.setFromOptions(); CHKERRQ(ierr);
      ierr = iceFactory.create(&ice); CHKERRQ(ierr);
    }

    // set options specific to this particular ice type:
    ierr = ice->setFromOptions(); CHKERRQ(ierr);
  }

  // Create the stress balance object:
  PetscScalar pseudo_plastic_q = config.get("pseudo_plastic_q"),
    pseudo_plastic_uthreshold = config.get("pseudo_plastic_uthreshold") / secpera,
    plastic_regularization = config.get("plastic_regularization") / secpera;

  bool do_pseudo_plastic_till = config.get_flag("do_pseudo_plastic_till");
  
  if (basal == NULL)
    basal = new IceBasalResistancePlasticLaw(plastic_regularization, do_pseudo_plastic_till, 
                                             pseudo_plastic_q, pseudo_plastic_uthreshold);

  if (EC == NULL) {
    EC = new EnthalpyConverter(config);
    if (getVerbosityLevel() > 3) {
      PetscViewer viewer;
      ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
      ierr = EC->viewConstants(viewer); CHKERRQ(ierr);
    }
  }

  // If both SIA and SSA are "on", the SIA and SSA velocities are always added
  // up (there is no switch saying "do the hybrid").
  if (stress_balance == NULL) {
    ShallowStressBalance *my_stress_balance;

    SSB_Modifier *modifier = new SIAFD(grid, *ice, *EC, config);

    if (expername == 'G' || expername == 'H') {
      my_stress_balance = new SIA_Sliding(grid, *basal, *ice, *EC, config);
    } else {
      my_stress_balance = new SSB_Trivial(grid, *basal, *ice, *EC, config);
    }
  
    // ~PISMStressBalance() will de-allocate my_stress_balance and modifier.
    stress_balance = new PISMStressBalance(grid, my_stress_balance,
                                           modifier, NULL, config);

    // Note that in PISM stress balance computations are diagnostic, i.e. do not
    // have a state that changes in time. This means that this call can be here
    // and not in model_state_setup() and we don't need to re-initialize after
    // the "diagnostic time step".
    ierr = stress_balance->init(variables); CHKERRQ(ierr);

    if (config.get_flag("include_bmr_in_continuity")) {
      ierr = stress_balance->set_basal_melt_rate(&vbmr); CHKERRQ(ierr);
    }
  }

  // see EISMINT II description; choose no ocean interaction, purely SIA, and E=1
  config.set_flag("is_dry_simulation", true);
  config.set_flag("use_ssa_velocity", false);

  // Make bedrock thermal material properties into ice properties.  Note that
  // zero thickness bedrock layer is the default, but we want the ice/rock
  // interface segment to have geothermal flux applied directly to ice without
  // jump in material properties at base.
  // (These must follow call to IceModel::init_physics(), where ice is allocated
  // as a pointer.)
  config.set("bedrock_thermal_density", ice->rho);
  config.set("bedrock_thermal_conductivity", ice->k);
  config.set("bedrock_thermal_specific_heat_capacity", ice->c_p);

  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceEISModel::init_couplers() {
  PetscErrorCode      ierr;

  ierr = IceModel::init_couplers(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
    "  setting surface mass balance and surface temperature variables from formulas...\n");
  CHKERRQ(ierr);

  // now fill in accum and surface temp
  ierr = artm.begin_access(); CHKERRQ(ierr);
  ierr = acab.begin_access(); CHKERRQ(ierr);

  PetscScalar cx = grid.Lx, cy = grid.Ly;
  if (expername == 'E') {  cx += 100.0e3;  cy += 100.0e3;  } // shift center
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // r is distance from center of grid; if E then center is shifted (above)
      const PetscScalar r = sqrt( PetscSqr(-cx + grid.dx*i)
                                  + PetscSqr(-cy + grid.dy*j) );
      // set accumulation from formula (7) in (Payne et al 2000)
      acab(i,j) = PetscMin(M_max, S_b * (R_el-r));
      // set surface temperature
      artm(i,j) = T_min + S_T * r;  // formula (8) in (Payne et al 2000)
    }
  }

  ierr = artm.end_access(); CHKERRQ(ierr);
  ierr = acab.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::generateTroughTopography() {
  PetscErrorCode  ierr;
  // computation based on code by Tony Payne, 6 March 1997:
  //    http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  
  const PetscScalar    b0 = 1000.0;  // plateau elevation
  const PetscScalar    L = 750.0e3;  // half-width of computational domain
  const PetscScalar    w = 200.0e3;  // trough width
  const PetscScalar    slope = b0/L;
  const PetscScalar    dx61 = (2*L) / 60; // = 25.0e3
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * grid.dx, ewd = j * grid.dy;
      if (    (nsd >= (27 - 1) * dx61) && (nsd <= (35 - 1) * dx61)
           && (ewd >= (31 - 1) * dx61) && (ewd <= (61 - 1) * dx61) ) {
        vbed(i,j) = 1000.0 - PetscMax(0.0, slope * (ewd - L) * cos(pi * (nsd - L) / w));
      } else {
        vbed(i,j) = 1000.0;
      }
    }
  }
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
               "trough bed topography stored by IceEISModel::generateTroughTopography()\n");
               CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceEISModel::generateMoundTopography() {
  PetscErrorCode  ierr;
  // computation based on code by Tony Payne, 6 March 1997:
  //    http://homepages.vub.ac.be/~phuybrec/eismint/topog2.f
  
  const PetscScalar    slope = 250.0;
  const PetscScalar    w = 150.0e3;  // mound width
  ierr = vbed.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar nsd = i * grid.dx, ewd = j * grid.dy;
      vbed(i,j) = PetscAbs(slope * sin(pi * ewd / w) + slope * cos(pi * nsd / w));
    }
  }
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,
    "mound bed topography stored by IceEISModel::generateTroughTopography()\n");
    CHKERRQ(ierr);
  return 0;
}


//! Only executed if NOT initialized from file (-i).
PetscErrorCode IceEISModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // initialize from EISMINT II formulas
  ierr = verbPrintf(2,grid.com, 
    "initializing variables from EISMINT II experiment %c formulas ... \n", 
    expername); CHKERRQ(ierr);

  ierr = vbed.set(0.0);
  if ((expername == 'I') || (expername == 'J')) {
    ierr = generateTroughTopography(); CHKERRQ(ierr);
  } 
  if ((expername == 'K') || (expername == 'L')) {
    ierr = generateMoundTopography(); CHKERRQ(ierr);
  } 
  // communicate b in any case; it will be horizontally-differentiated
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);

  ierr = vHmelt.set(0.0); CHKERRQ(ierr);
  ierr = vbmr.set(0.0); CHKERRQ(ierr);
  ierr = vGhf.set(0.042); CHKERRQ(ierr);  // EISMINT II value; J m-2 s-1

  ierr = vMask.set(MASK_GROUNDED); CHKERRQ(ierr);
  ierr = vuplift.set(0.0); CHKERRQ(ierr);  // no expers have uplift at start

  ierr = vtillphi.set(config.get("default_till_phi")); CHKERRQ(ierr);

  // if no -i file then starts with zero ice
  ierr = vh.set(0.0); CHKERRQ(ierr);
  ierr = vH.set(0.0); CHKERRQ(ierr);

  ierr = regrid(2); CHKERRQ(ierr);
  
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  // this IceModel bootstrap method should do right thing because of variable
  //   settings above and init of coupler above
  ierr = putTempAtDepth(); CHKERRQ(ierr);

  ierr = regrid(3); CHKERRQ(ierr);

  return 0;
}

