// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

//! \file pismo.cc A regional (outlet glacier) model form of PISM.
/*! \file pismo.cc 
The classes in this file modify basic PISM whole ice sheet modeling assumptions.
Normally in PISM the ice sheet occupies a continent which is surrounded by
ocean.  Or at least PISM assumes that the edge of the computational domain is in
a region with strong ablation that the ice will not cross.

Here, by contrast, we add a strip around the edge of the computational domain
(variable `no_model_mask` and option `-no_model_strip`).  Various
simplifications and boundary conditions are enforced in this script:
* the surface gradient computation is made trivial,
* the driving stress does not change during the run but instead comes from
the gradient of a saved surface elevation, and
* the base is made strong so that no sliding occurs.

Also options `-force_to_thickness_file` and variable `ftt_mask` play a role in isolating
the modeled outlet glacier.  But there is no code here for that purpose. 
Instead see the PSForceThickness surface model modifier class.
 */

//! \brief A version of the PISM core class (IceModel) which knows about the
//! no_model_mask and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, Config &c, Config &o)
     : IceModel(g,c,o) {};
protected:
  virtual void set_vars_from_options();
  virtual void bootstrap_2d(const std::string &filename);
  virtual void initFromFile(const std::string &filename);
  virtual void model_state_setup();
  virtual void createVecs();
  virtual void allocate_stressbalance();
  virtual void allocate_basal_yield_stress();
  virtual void massContExplicitStep();
  virtual void cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     StarStencil<Vector2> input_velocity,
                                     StarStencil<double> input_flux,
                                     StarStencil<double> &output_velocity,
                                     StarStencil<double> &output_flux);
  virtual void enthalpyAndDrainageStep(double* vertSacrCount,
                                                 double* liquifiedVol,
                                                 double* bulgeCount);
private:
  IceModelVec2Int no_model_mask;
  IceModelVec2S   usurfstore, thkstore;
  IceModelVec2S   bmr_stored;
  void  set_no_model_strip(double stripwidth);
};

//! \brief Set no_model_mask variable to have value 1 in strip of width 'strip'
//! m around edge of computational domain, and value 0 otherwise.
void IceRegionalModel::set_no_model_strip(double strip) {

  IceModelVec::AccessList list(no_model_mask);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (in_null_strip(grid, i, j, strip) == true) {
      no_model_mask(i, j) = 1;
    } else {
      no_model_mask(i, j) = 0;
    }
  }

  no_model_mask.metadata().set_string("pism_intent", "model_state");

  no_model_mask.update_ghosts();
}


void IceRegionalModel::createVecs() {

  IceModel::createVecs();

  verbPrintf(2, grid.com,
             "  creating IceRegionalModel vecs ...\n");

  // stencil width of 2 needed for surfaceGradientSIA() action
  no_model_mask.create(grid, "no_model_mask", WITH_GHOSTS, 2);
  no_model_mask.set_attrs("model_state", // ensures that it gets written at the end of the run
                          "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)",
                          "", ""); // no units and no standard name
  double NMMASK_NORMAL   = 0.0,
         NMMASK_ZERO_OUT = 1.0;
  std::vector<double> mask_values(2);
  mask_values[0] = NMMASK_NORMAL;
  mask_values[1] = NMMASK_ZERO_OUT;
  no_model_mask.metadata().set_doubles("flag_values", mask_values);
  no_model_mask.metadata().set_string("flag_meanings", "normal special_treatment");
  no_model_mask.set_time_independent(true);
  no_model_mask.set(NMMASK_NORMAL);
  grid.variables().add(no_model_mask);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  usurfstore.create(grid, "usurfstore", WITH_GHOSTS, 2);
  usurfstore.set_attrs("model_state", // ensures that it gets written at the end of the run
                       "saved surface elevation for use to keep surface gradient constant in no_model strip",
                       "m",
                       ""); //  no standard name
  grid.variables().add(usurfstore);

  // stencil width of 1 needed for differentiation
  thkstore.create(grid, "thkstore", WITH_GHOSTS, 1);
  thkstore.set_attrs("model_state", // ensures that it gets written at the end of the run
                     "saved ice thickness for use to keep driving stress constant in no_model strip",
                     "m",
                     ""); //  no standard name
  grid.variables().add(thkstore);

  // Note that the name of this variable (bmr_stored) does not matter: it is
  // *never* read or written. We make a copy of bmelt instead.
  bmr_stored.create(grid, "bmr_stored", WITH_GHOSTS, 2);
  bmr_stored.set_attrs("internal",
                       "time-independent basal melt rate in the no-model-strip",
                       "m s-1", "");

  if (config.get_flag("ssa_dirichlet_bc")) {
    // remove the bcflag variable from the dictionary
    grid.variables().remove("bcflag");

    grid.variables().add(no_model_mask, "bcflag");
  }
}

void IceRegionalModel::model_state_setup() {

  IceModel::model_state_setup();

  // Now save the basal melt rate at the beginning of the run.
  bmr_stored.copy_from(basal_melt_rate);

  bool zgwnm = options::Bool("-zero_grad_where_no_model",
                             "set zero surface gradient in no model strip");
  if (zgwnm) {
    thkstore.set(0.0);
    usurfstore.set(0.0);
  }

  options::Real strip_km("-no_model_strip", 
                        "width in km of strip near boundary in which modeling is turned off",
                        0.0);

  if (strip_km.is_set()) {
    verbPrintf(2, grid.com,
               "* Option -no_model_strip read... setting boundary strip width to %.2f km\n",
               strip_km.value());
    set_no_model_strip(grid.convert(strip_km, "km", "m"));
  }
}

void IceRegionalModel::allocate_stressbalance() {

  if (stress_balance != NULL) {
    return;
  }

  std::string model = config.get_string("stress_balance_model");

  ShallowStressBalance *sliding = NULL;
  if (model == "none" || model == "sia") {
    sliding = new ZeroSliding(grid, *EC);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(grid, *EC);
  } else if (model == "ssa" || model == "ssa+sia") {
    sliding = new SSAFD_Regional(grid, *EC);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(grid, *EC);
  } else if (model == "prescribed_sliding+sia" || "ssa+sia") {
    modifier = new SIAFD_Regional(grid, *EC);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  // ~StressBalance() will de-allocate sliding and modifier.
  stress_balance = new StressBalance(grid, sliding, modifier);
}


void IceRegionalModel::allocate_basal_yield_stress() {

  if (basal_yield_stress_model != NULL) {
    return;
  }

  std::string model = config.get_string("stress_balance_model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = config.get_string("yield_stress_model");

    if (yield_stress_model == "constant") {
      basal_yield_stress_model = new ConstantYieldStress(grid);
    } else if (yield_stress_model == "mohr_coulomb") {
      basal_yield_stress_model = new RegionalDefaultYieldStress(grid, subglacial_hydrology);
    } else {
      throw RuntimeError::formatted("yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }
  }
}


void IceRegionalModel::bootstrap_2d(const std::string &filename) {

  IceModel::bootstrap_2d(filename);

  usurfstore.regrid(filename, OPTIONAL, 0.0);
  thkstore.regrid(filename, OPTIONAL, 0.0);
}


void IceRegionalModel::initFromFile(const std::string &filename) {
  PIO nc(grid, "guess_mode");

  bool no_model_strip_set = options::Bool("-no_model_strip", "No-model strip, in km");

  if (no_model_strip_set) {
    no_model_mask.metadata().set_string("pism_intent", "internal");
  }

  verbPrintf(2, grid.com,
             "* Initializing IceRegionalModel from NetCDF file '%s'...\n",
             filename.c_str());

  // Allow re-starting from a file that does not contain u_ssa_bc and v_ssa_bc.
  // The user is probably using -regrid_file to bring in SSA B.C. data.
  if (config.get_flag("ssa_dirichlet_bc")) {
    bool u_ssa_exists, v_ssa_exists;

    nc.open(filename, PISM_READONLY);
    u_ssa_exists = nc.inq_var("u_ssa_bc");
    v_ssa_exists = nc.inq_var("v_ssa_bc");
    nc.close();

    if (! (u_ssa_exists && v_ssa_exists)) {
      vBCvel.metadata().set_string("pism_intent", "internal");
      verbPrintf(2, grid.com,
                 "PISM WARNING: u_ssa_bc and/or v_ssa_bc not found in %s. Setting them to zero.\n"
                 "              This may be overridden by the -regrid_file option.\n",
                 filename.c_str());

      vBCvel.set(0.0);
    }
  }

  bool zgwnm = options::Bool("-zero_grad_where_no_model",
                             "zero surface gradient in no model strip");
  if (zgwnm) {
    thkstore.metadata().set_string("pism_intent", "internal");
    usurfstore.metadata().set_string("pism_intent", "internal");
  }

  IceModel::initFromFile(filename);

  if (config.get_flag("ssa_dirichlet_bc")) {
      vBCvel.metadata().set_string("pism_intent", "model_state");
  }

  if (zgwnm) {
    thkstore.metadata().set_string("pism_intent", "model_state");
    usurfstore.metadata().set_string("pism_intent", "model_state");
  }
}


void IceRegionalModel::set_vars_from_options() {

  // base class reads the -boot_file option and does the bootstrapping:
  IceModel::set_vars_from_options();

  bool nmstripSet = options::Bool("-no_model_strip", 
                                 "width in km of strip near boundary in which modeling is turned off");

  if (not nmstripSet) {
    throw RuntimeError("option '-no_model_strip X' (X in km) is REQUIRED if '-i' is not used.\n"
                       "pismo has no well-defined semantics without it!");
  }

  if (config.get_flag("do_cold_ice_methods")) {
    throw RuntimeError("pismo does not support the 'cold' mode.");
  }
}

void IceRegionalModel::massContExplicitStep() {

  // This ensures that no_model_mask is available in
  // IceRegionalModel::cell_interface_fluxes() below.
  IceModelVec::AccessList list(no_model_mask);

  IceModel::massContExplicitStep();
}

void IceRegionalModel::cell_interface_fluxes(bool dirichlet_bc,
                                             int i, int j,
                                             StarStencil<Vector2> input_velocity,
                                             StarStencil<double> input_flux,
                                             StarStencil<double> &output_velocity,
                                             StarStencil<double> &output_flux) {

  IceModel::cell_interface_fluxes(dirichlet_bc, i, j,
                                  input_velocity,
                                  input_flux,
                                  output_velocity,
                                  output_flux);

  StarStencil<int> nmm = no_model_mask.int_star(i,j);
  Direction dirs[4] = {North, East, South, West};

  for (int n = 0; n < 4; ++n) {
    Direction direction = dirs[n];

      if ((nmm.ij == 1) || (nmm.ij == 0 && nmm[direction] == 1)) {
      output_velocity[direction] = 0.0;
      output_flux[direction] = 0.0;
    }
  }
  //
}

void IceRegionalModel::enthalpyAndDrainageStep(double* vertSacrCount, double* liquifiedVol,
                                                         double* bulgeCount) {
  double *new_enthalpy, *old_enthalpy;

  IceModel::enthalpyAndDrainageStep(vertSacrCount, liquifiedVol, bulgeCount);

  // note that the call above sets vWork3d; ghosts are comminucated later (in
  // IceModel::energyStep()).
  IceModelVec::AccessList list;
  list.add(no_model_mask);
  list.add(vWork3d);
  list.add(Enth3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (no_model_mask(i, j) < 0.5) {
      continue;
    }

    vWork3d.getInternalColumn(i, j, &new_enthalpy);
    Enth3.getInternalColumn(i, j, &old_enthalpy);

    for (unsigned int k = 0; k < grid.Mz(); ++k) {
      new_enthalpy[k] = old_enthalpy[k];
    }
  }

  // set basal_melt_rate; ghosts are comminucated later (in IceModel::energyStep()).
  list.add(basal_melt_rate);
  list.add(bmr_stored);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (no_model_mask(i, j) < 0.5) {
      continue;
    }

    basal_melt_rate(i, j) = bmr_stored(i, j);
  }
}


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  MPI_Comm com = MPI_COMM_WORLD;

  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();

    verbPrintf(2,com, "PISMO %s (regional outlet-glacier run mode)\n",
               PISM_Revision);

    if (options::Bool("-version", "stop after printing print PISM version")) {
      return 0;
    }

    bool iset = options::Bool("-i", "input file name");
    bool bfset = options::Bool("-boot_file", "bootstrapping file name");
    std::string usage =
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
    if ((not iset) && (not bfset)) {
      ierr = PetscPrintf(com,
                         "\nPISM ERROR: one of options -i,-boot_file is required\n\n");
      PISM_PETSC_CHK(ierr, "PetscPrintf");
      show_usage(com, "pismo", usage);
      return 0;
    } else {
      std::vector<std::string> required;
      required.clear();

      bool done = show_usage_check_req_opts(com, "pismo", required, usage);
      if (done) {
        return 0;
      }
    }

    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    init_config(com, config, overrides, true);

    // initialize the ice dynamics model
    IceGrid g(com, config);
    IceRegionalModel m(g, config, overrides);
    m.setExecName("pismo");

    m.init();

    m.run();

    verbPrintf(2,com, "... done with run\n");

    // provide a default output file name if no -o option is given.
    m.writeFiles("unnamed_regional.nc");
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}

