// Copyright (C) 2010, 2011, 2012 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

static char help[] =
  "Ice sheet driver for PISM regional (outlet glacier) simulations, initialized\n"
  "from data.\n";

#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"

#include "regional.hh"

#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"
#include "PISMStressBalance.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "PISMConstantYieldStress.hh"
#include "PIO.hh"
#include "pism_options.hh"

//! \file pismo.cc A regional (outlet glacier) model form of PISM.
/*! \file pismo.cc 
The classes in this file modify basic PISM whole ice sheet modeling assumptions.
Normally in PISM the ice sheet occupies a continent which is surrounded by
ocean.  Or at least PISM assumes that the edge of the computational domain is in
a region with strong ablation that the ice will not cross.

Here, by contrast, we add a strip around the edge of the computational domain
(variable \c no_model_mask and option \c -no_model_strip).  Various
simplifications and boundary conditions are enforced in this script:
* the surface gradient computation is made trivial,
* the driving stress does not change during the run but instead comes from
the gradient of a saved surface elevation, and
* the base is made strong so that no sliding occurs.

Also options \c -force_to_thk and variable \c ftt_mask play a role in isolating
the modeled outlet glacier.  But there is no code here for that purpose. 
Instead see the PSForceThickness surface model modifier class.
 */

//! \brief A version of the PISM core class (IceModel) which knows about the
//! no_model_mask and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o)
     : IceModel(g,c,o) {};
protected:
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode bootstrap_2d(string filename);
  virtual PetscErrorCode initFromFile(string filename);
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_basal_yield_stress();
  virtual PetscErrorCode massContExplicitStep();
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     planeStar<PISMVector2> input_velocity,
                                     planeStar<PetscScalar> input_flux,
                                     planeStar<PetscScalar> &output_velocity,
                                     planeStar<PetscScalar> &output_flux);
  virtual PetscErrorCode enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                                                 PetscScalar* bulgeCount);
  virtual PetscErrorCode setFromOptions();
private:
  IceModelVec2Int no_model_mask;
  IceModelVec2S   usurfstore, thkstore;
  IceModelVec2S   bmr_stored;
  PetscErrorCode  set_no_model_strip(PetscReal stripwidth);
};

//! \brief 
PetscErrorCode IceRegionalModel::setFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);

  ierr = config.flag_from_option("ssa_dirichlet_bc", "ssa_dirichlet_bc"); CHKERRQ(ierr);

  return 0;
}

//! \brief Set no_model_mask variable to have value 1 in strip of width 'strip'
//! m around edge of computational domain, and value 0 otherwise.
PetscErrorCode IceRegionalModel::set_no_model_strip(PetscReal strip) {
  PetscErrorCode ierr;

    ierr = no_model_mask.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (grid.x[i] <= grid.x[0]+strip || grid.x[i] >= grid.x[grid.Mx-1]-strip) {
          no_model_mask(i, j) = 1;
        } else if (grid.y[j] <= grid.y[0]+strip || grid.y[j] >= grid.y[grid.My-1]-strip) {
          no_model_mask(i, j) = 1;
        } else {
          no_model_mask(i, j) = 0;
        }
      }
    }
    ierr = no_model_mask.end_access(); CHKERRQ(ierr);

    ierr = no_model_mask.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);

    ierr = no_model_mask.beginGhostComm(); CHKERRQ(ierr);
    ierr = no_model_mask.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "  creating IceRegionalModel vecs ...\n"); CHKERRQ(ierr);

  // stencil width of 2 needed for surfaceGradientSIA() action
  ierr = no_model_mask.create(grid, "no_model_mask", true, 2); CHKERRQ(ierr);
  ierr = no_model_mask.set_attrs("model_state", // ensures that it gets written at the end of the run
    "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)",
    "", ""); CHKERRQ(ierr); // no units and no standard name
  double NMMASK_NORMAL   = 0.0,
         NMMASK_ZERO_OUT = 1.0;
  vector<double> mask_values(2);
  mask_values[0] = NMMASK_NORMAL;
  mask_values[1] = NMMASK_ZERO_OUT;
  ierr = no_model_mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = no_model_mask.set_attr("flag_meanings", "normal special_treatment"); CHKERRQ(ierr);
  no_model_mask.output_data_type = PISM_BYTE;
  no_model_mask.time_independent = true;
  ierr = no_model_mask.set(NMMASK_NORMAL); CHKERRQ(ierr);
  ierr = variables.add(no_model_mask); CHKERRQ(ierr);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  ierr = usurfstore.create(grid, "usurfstore", true, 2); CHKERRQ(ierr);
  ierr = usurfstore.set_attrs(
    "model_state", // ensures that it gets written at the end of the run
    "saved surface elevation for use to keep surface gradient constant in no_model strip",
    "m",
    ""); CHKERRQ(ierr); //  no standard name
  ierr = variables.add(usurfstore); CHKERRQ(ierr);

  // stencil width of 1 needed for differentiation
  ierr = thkstore.create(grid, "thkstore", true, 1); CHKERRQ(ierr);
  ierr = thkstore.set_attrs(
    "model_state", // ensures that it gets written at the end of the run
    "saved ice thickness for use to keep driving stress constant in no_model strip",
    "m",
    ""); CHKERRQ(ierr); //  no standard name
  ierr = variables.add(thkstore); CHKERRQ(ierr);

  // Note that the name of this variable (bmr_stored) does not matter: it is
  // *never* read or written. We make a copy of bmelt instead.
  ierr = bmr_stored.create(grid, "bmr_stored", true, 2); CHKERRQ(ierr);
  ierr = bmr_stored.set_attrs("internal",
                              "time-independent basal melt rate in the no-model-strip",
                              "m s-1", ""); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
    // remove the bcflag variable from the dictionary
    variables.remove("bcflag");

    ierr = variables.add(no_model_mask, "bcflag"); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceRegionalModel::model_state_setup() {
  PetscErrorCode ierr;

  ierr = IceModel::model_state_setup(); CHKERRQ(ierr);

  // Now save the basal melt rate at the beginning of the run.
  ierr = bmr_stored.copy_from(vbmr); CHKERRQ(ierr);

  bool zgwnm;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (zgwnm) {
    ierr = thkstore.set(0.0); CHKERRQ(ierr);
    ierr = usurfstore.set(0.0); CHKERRQ(ierr);
  }

  bool nmstripSet;
  PetscReal stripkm = 0.0;
  ierr = PISMOptionsReal("-no_model_strip", 
                         "width in km of strip near boundary in which modeling is turned off",
			 stripkm, nmstripSet);

  if (nmstripSet) {
    ierr = verbPrintf(2, grid.com,
                      "* Option -no_model_strip read... setting boundary strip width to %.2f km\n",
                      stripkm); CHKERRQ(ierr);
    ierr = set_no_model_strip(convert(stripkm, "km", "m")); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceRegionalModel::allocate_stressbalance() {
  PetscErrorCode ierr;
  
  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_sia = config.get_flag("do_sia");

  ShallowStressBalance *my_stress_balance;
  SSB_Modifier *modifier;
  if (do_sia) {
    modifier = new SIAFD_Regional(grid, *EC, config);
  } else {
    modifier = new SSBM_Trivial(grid, *EC, config);
  }

  if (use_ssa_velocity) {
    my_stress_balance = new SSAFD_Regional(grid, *basal, *EC, config);
  } else {
    my_stress_balance = new SSB_Trivial(grid, *basal, *EC, config);
  }
  
  // ~PISMStressBalance() will de-allocate my_stress_balance and modifier.
  stress_balance = new PISMStressBalance(grid, my_stress_balance,
                                         modifier, ocean, config);

  // Note that in PISM stress balance computations are diagnostic, i.e. do not
  // have a state that changes in time. This means that this call can be here
  // and not in model_state_setup() and we don't need to re-initialize after
  // the "diagnostic time step".
  ierr = stress_balance->init(variables); CHKERRQ(ierr);

  if (config.get_flag("include_bmr_in_continuity")) {
    ierr = stress_balance->set_basal_melt_rate(&vbmr); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceRegionalModel::allocate_basal_yield_stress() {
  PetscErrorCode ierr;

  if (basal_yield_stress != NULL)
    return 0;

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_blatter = config.get_flag("do_blatter");

  if (use_ssa_velocity || do_blatter) {
    bool hold_tauc;
    ierr = PISMOptionsIsSet("-hold_tauc", hold_tauc); CHKERRQ(ierr);
    
    if (hold_tauc) {
      basal_yield_stress = new PISMConstantYieldStress(grid, config);
    } else {
      basal_yield_stress = new PISMRegionalDefaultYieldStress(grid, config);
    }
  }
  
  return 0;
}


PetscErrorCode IceRegionalModel::bootstrap_2d(string filename) {
  PetscErrorCode ierr;

  ierr = IceModel::bootstrap_2d(filename); CHKERRQ(ierr);

  ierr = usurfstore.regrid(filename, 0.0); CHKERRQ(ierr);
  ierr =   thkstore.regrid(filename, 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceRegionalModel::initFromFile(string filename) {
  PetscErrorCode  ierr;
  PIO nc(grid, "guess_format");

  bool no_model_strip_set;
  ierr = PISMOptionsIsSet("-no_model_strip", "No-model strip, in km",
                          no_model_strip_set); CHKERRQ(ierr);

  if (no_model_strip_set) {
    ierr = no_model_mask.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
                    "* Initializing IceRegionalModel from NetCDF file '%s'...\n",
                    filename.c_str()); CHKERRQ(ierr);

  // Allow re-starting from a file that does not contain u_ssa_bc and v_ssa_bc.
  // The user is probably using -regrid_file to bring in SSA B.C. data.
  if (config.get_flag("ssa_dirichlet_bc")) {
    bool u_ssa_exists, v_ssa_exists;

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("u_ssa_bc", u_ssa_exists); CHKERRQ(ierr);
    ierr = nc.inq_var("v_ssa_bc", v_ssa_exists); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    if (! (u_ssa_exists && v_ssa_exists)) {
      ierr = vBCvel.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
      ierr = verbPrintf(2, grid.com,
                        "PISM WARNING: u_ssa_bc and/or v_ssa_bc not found in %s. Setting them to zero.\n"
                        "              This may be overridden by the -regrid_file option.\n",
                        filename.c_str()); CHKERRQ(ierr);

      ierr = vBCvel.set(0.0); CHKERRQ(ierr);
    }
  }

  bool zgwnm;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (zgwnm) {
    ierr = thkstore.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
    ierr = usurfstore.set_attr("pism_intent", "internal"); CHKERRQ(ierr);
  }

  ierr = IceModel::initFromFile(filename); CHKERRQ(ierr);

  if (config.get_flag("ssa_dirichlet_bc")) {
      ierr = vBCvel.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  }

  if (zgwnm) {
    ierr = thkstore.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
    ierr = usurfstore.set_attr("pism_intent", "model_state"); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceRegionalModel::set_vars_from_options() {
  PetscErrorCode ierr;
  bool nmstripSet;

  // base class reads the -boot_file option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  ierr = PISMOptionsIsSet("-no_model_strip", 
                          "width in km of strip near boundary in which modeling is turned off",
                          nmstripSet);
  if (!nmstripSet) {
    ierr = PetscPrintf(grid.com,
      "PISMO ERROR: option '-no_model_strip X' (X in km) is REQUIRED if '-i' is not used.\n"
      "   pismo has no well-defined semantics without it!  ENDING ...\n\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  if (config.get_flag("do_cold_ice_methods")) {
    PetscPrintf(grid.com, "PISM ERROR: pismo does not support the 'cold' mode.\n");
    PISMEnd();
  }

  return 0;
}

PetscErrorCode IceRegionalModel::massContExplicitStep() {
  PetscErrorCode ierr;

  // This ensures that no_model_mask is available in
  // IceRegionalModel::cell_interface_fluxes() below.
  ierr = no_model_mask.begin_access(); CHKERRQ(ierr);

  ierr = IceModel::massContExplicitStep(); CHKERRQ(ierr);

  ierr = no_model_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

void IceRegionalModel::cell_interface_fluxes(bool dirichlet_bc,
                                             int i, int j,
                                             planeStar<PISMVector2> input_velocity,
                                             planeStar<PetscScalar> input_flux,
                                             planeStar<PetscScalar> &output_velocity,
                                             planeStar<PetscScalar> &output_flux) {

  IceModel::cell_interface_fluxes(dirichlet_bc, i, j,
                                  input_velocity,
                                  input_flux,
                                  output_velocity,
                                  output_flux);

  planeStar<int> nmm = no_model_mask.int_star(i,j);
  PISM_Direction dirs[4] = {North, East, South, West};

  for (int n = 0; n < 4; ++n) {
    PISM_Direction direction = dirs[n];

      if ((nmm.ij == 1) ||
          (nmm.ij == 0 && nmm[direction] == 1)) {
      output_velocity[direction] = 0.0;
      output_flux[direction] = 0.0;
    }
  }
  //
}

PetscErrorCode IceRegionalModel::enthalpyAndDrainageStep(PetscScalar* vertSacrCount, PetscScalar* liquifiedVol,
                                                         PetscScalar* bulgeCount) {
  PetscErrorCode ierr;
  PetscScalar *new_enthalpy, *old_enthalpy;

  ierr = IceModel::enthalpyAndDrainageStep(vertSacrCount, liquifiedVol, bulgeCount); CHKERRQ(ierr);

  // note that the call above sets vWork3d; ghosts are comminucated later (in
  // IceModel::energyStep()).
  ierr = no_model_mask.begin_access(); CHKERRQ(ierr);

  ierr = vWork3d.begin_access(); CHKERRQ(ierr);
  ierr = Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (no_model_mask(i, j) < 0.5)
        continue;

      ierr = vWork3d.getInternalColumn(i, j, &new_enthalpy); CHKERRQ(ierr);
      ierr = Enth3.getInternalColumn(i, j, &old_enthalpy); CHKERRQ(ierr);

      for (int k = 0; k < grid.Mz; ++k)
        new_enthalpy[k] = old_enthalpy[k];
    }
  }
  ierr = Enth3.end_access(); CHKERRQ(ierr);
  ierr = vWork3d.end_access(); CHKERRQ(ierr);

  // set vbmr; ghosts are comminucated later (in IceModel::energyStep()).
  ierr = vbmr.begin_access(); CHKERRQ(ierr);
  ierr = bmr_stored.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (no_model_mask(i, j) < 0.5)
        continue;

      vbmr(i, j) = bmr_stored(i, j);
    }
  }
  ierr = bmr_stored.end_access(); CHKERRQ(ierr);
  ierr = vbmr.end_access(); CHKERRQ(ierr);

  ierr = no_model_mask.end_access(); CHKERRQ(ierr);

  return 0;
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  MPI_Comm    com = PETSC_COMM_WORLD;
  PetscMPIInt rank, size;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "PISMO %s (regional outlet-glacier run mode)\n",
		      PISM_Revision); CHKERRQ(ierr);
    ierr = stop_on_version_option(); CHKERRQ(ierr);

    bool iset, bfset;
    ierr = PISMOptionsIsSet("-i", iset); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", bfset); CHKERRQ(ierr);
    string usage =
      "  pismo {-i IN.nc|-boot_file IN.nc} [-no_model_strip X] [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "  -no_model_strip X (re-)set width of no-model strip along edge of\n"
      "              computational domain to X km\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((!iset) && (!bfset)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismo", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;
      required.clear();
      ierr = show_usage_check_req_opts(com, "pismo", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides, true); CHKERRQ(ierr);

    // initialize the ice dynamics model
    IceGrid g(com, rank, size, config);
    IceRegionalModel m(g, config, overrides);
    ierr = m.setExecName("pismo"); CHKERRQ(ierr);

    // initialize boundary models
    // factories allow runtime choice
    PAFactory pa(g, config);
    PSFactory ps(g, config);
    POFactory po(g, config);
    // now read options and choose
    PISMAtmosphereModel *atmosphere;
    PISMSurfaceModel    *surface;
    PISMOceanModel      *ocean;
    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    surface->attach_atmosphere_model(atmosphere); // IceModel m does not see atmosphere
    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed_regional.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

