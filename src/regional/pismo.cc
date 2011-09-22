// Copyright (C) 2010, 2011 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

#include "PCFactory.hh"
#include "PISMAtmosphere.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMStressBalance.hh"
#include "SIAFD.hh"
#include "SSAFD.hh"

//! \file pismo.cc A regional (outlet glacier) model form of PISM.
/*! \file pismo.cc 
The classes in this file modify basic PISM whole ice sheet modeling
assumptions.  Especially that the ice sheet occupies a
continent which is surrounded by ocean.  (Or that the edge of the computational
domain is in a region with strong ablation that the ice will not cross.)

Various simplifications and boundary conditions are enforced in a strip around
the edge of the computational domain (variable \c no_model_mask and option
\c -no_model_strip):
* the surface gradient computation is made trivial
* the driving stress changes in the same way
* the base is made strong so no sliding occurs.

Also options \c -force_to_thk and variable \c ftt_mask play a role in isolating
the modeled outlet glacier.  See the PSForceThickness surface model modifier 
class.
 */

//! \brief A version of the SIA stress balance with tweaks for outlet glacier
//! simulations.
class SIAFD_Regional : public SIAFD
{
public:
  SIAFD_Regional(IceGrid &g, IceFlowLaw &i, EnthalpyConverter &e, const NCConfigVariable &c)
    : SIAFD(g, i, e, c) {}
  virtual ~SIAFD_Regional() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
protected:
  IceModelVec2Int *no_model_mask;
  IceModelVec2S   *usurfstore;   
};

PetscErrorCode SIAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SIAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SIA solver...\n"); CHKERRQ(ierr);

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(1, "no_model_mask is not available");

  usurfstore = dynamic_cast<IceModelVec2S*>(vars.get("usurfstore"));
  if (usurfstore == NULL) SETERRQ(1, "usurfstore is not available");

  return 0;
}

PetscErrorCode SIAFD_Regional::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  ierr = SIAFD::compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  IceModelVec2Int nmm = *no_model_mask;
  IceModelVec2S hst = *usurfstore; // convenience

  const int Mx = grid.Mx, My = grid.My;
  const PetscScalar dx = grid.dx, dy = grid.dy;  // convenience

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = nmm.begin_access(); CHKERRQ(ierr);
  ierr = hst.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 1;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

      // x-component, i-offset
      if (nmm(i, j) > 0.5 || nmm(i + 1, j) > 0.5) {

        if (i < 0 || i + 1 > Mx - 1)
          h_x(i, j, 0) = 0.0;
        else
          h_x(i, j, 0) = (hst(i + 1, j) - hst(i, j)) / dx;
      }

      // x-component, j-offset
      if (nmm(i - 1, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
          nmm(i - 1, j)     > 0.5 || nmm(i + 1, j)     > 0.5) {

        if (i - 1 < 0 || j + 1 > My - 1 || i + 1 > Mx - 1)
          h_x(i, j, 1) = 0.0;
        else
          h_x(i, j, 1) = ( + hst(i + 1, j + 1) + hst(i + 1, j)
                           - hst(i - 1, j + 1) - hst(i - 1, j) ) / (4.0 * dx);

      }

      // y-component, i-offset
      if (nmm(i, j + 1) > 0.5 || nmm(i + 1, j + 1) > 0.5 ||
          nmm(i, j - 1) > 0.5 || nmm(i + 1, j - 1) > 0.5) {
        if (i < 0 || j + 1 > My - 1 || i + 1 > Mx - 1 || j - 1 < 0)
          h_y(i, j, 0) = 0.0;
        else
          h_y(i, j, 0) = ( + hst(i + 1, j + 1) + hst(i, j + 1)
                           - hst(i + 1, j - 1) - hst(i, j - 1) ) / (4.0 * dy);
      }

      // y-component, j-offset
      if (nmm(i, j) > 0.5 || nmm(i, j + 1) > 0.5) {
        
        if (j < 0 || j + 1 > My - 1)
          h_y(i, j, 1) = 0.0;
        else
          h_y(i, j, 1) = (hst(i, j + 1) - hst(i, j)) / dy;
      }

    }
  }
  ierr = nmm.end_access(); CHKERRQ(ierr);
  ierr = hst.end_access(); CHKERRQ(ierr);
  ierr = h_y.end_access(); CHKERRQ(ierr);
  ierr = h_x.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief A version of the SSA stress balance with tweaks for outlet glacier
//! simulations.
class SSAFD_Regional : public SSAFD
{
public:
  SSAFD_Regional(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i, EnthalpyConverter &e,
                 const NCConfigVariable &c)
    : SSAFD(g, b, i, e, c) {}
  virtual ~SSAFD_Regional() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode compute_driving_stress(IceModelVec2V &taud);
protected:
  IceModelVec2Int *no_model_mask;    
  IceModelVec2S   *usurfstore, *thkstore;
};

PetscErrorCode SSAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SSA solver...\n"); CHKERRQ(ierr);

  if (config.get_flag("dirichlet_bc")) {
    ierr = verbPrintf(2,grid.com,"  using stored SSA velocities as Dirichlet B.C. in the no_model_strip...\n"); 
    CHKERRQ(ierr);
  }

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(1, "no_model_mask is not available");

  usurfstore = dynamic_cast<IceModelVec2S*>(vars.get("usurfstore"));
  if (usurfstore == NULL) SETERRQ(1, "usurfstore is not available");

  thkstore = dynamic_cast<IceModelVec2S*>(vars.get("thkstore"));
  if (thkstore == NULL) SETERRQ(1, "thkstore is not available");

  return 0;
}

PetscErrorCode SSAFD_Regional::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  ierr = SSAFD::compute_driving_stress(result); CHKERRQ(ierr);

  const PetscReal standard_gravity = config.get("standard_gravity");
  IceModelVec2Int nmm = *no_model_mask;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = nmm.begin_access(); CHKERRQ(ierr);
  ierr = usurfstore->begin_access(); CHKERRQ(ierr);
  ierr = thkstore->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar pressure = ice.rho * standard_gravity * (*thkstore)(i,j);
      if (pressure <= 0) pressure = 0;

      if (nmm(i, j) > 0.5 || nmm(i - 1, j) > 0.5 || nmm(i + 1, j) > 0.5) {
        if (i - 1 < 0 || i + 1 > grid.Mx - 1)
          result(i, j).u = 0;
        else
          result(i, j).u = - pressure * usurfstore->diff_x(i,j);
      }

      if (nmm(i, j) > 0.5 || nmm(i, j - 1) > 0.5 || nmm(i, j + 1) > 0.5) {
        if (j - 1 < 0 || j + 1 > grid.My - 1)
          result(i, j).v = 0;
        else
          result(i, j).v = - pressure * usurfstore->diff_y(i,j);
      }

    }
  }
  ierr = usurfstore->end_access(); CHKERRQ(ierr);
  ierr = thkstore->end_access(); CHKERRQ(ierr);
  ierr = nmm.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


class PISMRegionalDefaultYieldStress : public PISMDefaultYieldStress
{
public:
  PISMRegionalDefaultYieldStress(IceGrid &g, const NCConfigVariable &conf)
    : PISMDefaultYieldStress(g, conf) {}
  virtual ~PISMRegionalDefaultYieldStress() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  IceModelVec2Int *no_model_mask;
};


PetscErrorCode PISMRegionalDefaultYieldStress::init(PISMVars &vars) {
  PetscErrorCode ierr;
  PetscInt v = getVerbosityLevel(); // turn off second, redundant init message
  ierr = setVerbosityLevel(1); CHKERRQ(ierr);
  ierr = PISMDefaultYieldStress::init(vars); CHKERRQ(ierr);
  ierr = setVerbosityLevel(v); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
    "  using the regional version with strong till in no_model_mask==1 area ...\n");
    CHKERRQ(ierr);
  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(1, "no_model_mask is not available");
  return 0;
}


PetscErrorCode PISMRegionalDefaultYieldStress::basal_material_yield_stress(IceModelVec2S &result) {
  PetscErrorCode ierr;
  
  // do whatever you normally do
  ierr = PISMDefaultYieldStress::basal_material_yield_stress(result); CHKERRQ(ierr);

  // now set result=tauc to a big value in no_model_strip
  ierr = no_model_mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*no_model_mask)(i,j) > 0.5) {
        result(i,j) = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar
      }
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = no_model_mask->end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief A version of the PISM core class (IceModel) which knows about the
//! no_model_mask and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o)
     : IceModel(g,c,o) {};
protected:
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode bootstrap_2d(const char *filename);
  virtual PetscErrorCode initFromFile(const char *filename);
  virtual PetscErrorCode model_state_setup();
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode allocate_stressbalance();
  virtual PetscErrorCode allocate_basal_yield_stress();
  virtual PetscErrorCode massContExplicitStep();
  virtual PetscErrorCode cell_interface_velocities(bool do_part_grid,
                                                   int i, int j,
                                                   planeStar<PetscScalar> &vel_output);
  virtual PetscErrorCode cell_interface_diffusive_flux(IceModelVec2Stag &Qstag, int i, int j,
                                                       planeStar<PetscScalar> &Q_output);
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

  ierr = config.flag_from_option("ssa_dirichlet_bc", "dirichlet_bc"); CHKERRQ(ierr);

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

    ierr = no_model_mask.beginGhostComm(); CHKERRQ(ierr);
    ierr = no_model_mask.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "  creating IceRegionalModel vecs ...\n"); CHKERRQ(ierr);

  // Most PISM variables are initialized in IceModel::model_state_setup() and
  // *not* in IceModel::createVecs().
  // 
  // We can get away with initializing no_model_mask, thkstore, usurfstore,
  // bmr_stored and enthalpy_stored here because all these variables are
  // constant in time.

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
  no_model_mask.output_data_type = NC_BYTE;
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
  ierr = usurfstore.set(0.0); CHKERRQ(ierr);
  ierr = variables.add(usurfstore); CHKERRQ(ierr);

  // stencil width of 1 needed for differentiation
  ierr = thkstore.create(grid, "thkstore", true, 1); CHKERRQ(ierr);
  ierr = thkstore.set_attrs(
    "model_state", // ensures that it gets written at the end of the run
    "saved ice thickness for use to keep driving stress constant in no_model strip",
    "m",
    ""); CHKERRQ(ierr); //  no standard name
  ierr = thkstore.set(0.0); CHKERRQ(ierr);
  ierr = variables.add(thkstore); CHKERRQ(ierr);

  // Note that the name of this variable (bmr_stored) does not matter: it is
  // *never* read or written. We make a copy of bmelt instead.
  ierr = bmr_stored.create(grid, "bmr_stored", true, 2); CHKERRQ(ierr);
  ierr = bmr_stored.set_attrs("internal",
                              "time-independent basal melt rate in the no-model-strip",
                              "m s-1", ""); CHKERRQ(ierr);

  if (config.get_flag("dirichlet_bc")) {
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

  if (config.get_flag("dirichlet_bc")) {
    double fill = convert(2e6, "m/year", "m/second");

    ierr = vBCMask.begin_access(); CHKERRQ(ierr);
    ierr = no_model_mask.begin_access(); CHKERRQ(ierr);
    ierr = vBCvel.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

        if (no_model_mask(i, j) < 0.5) {
          vBCvel(i, j).u = fill;
          vBCvel(i, j).v = fill;
          vBCMask(i, j)  = 0;
        } else {
          vBCMask(i, j) = 1;
        }

      }
    }
    ierr = vBCMask.end_access(); CHKERRQ(ierr);
    ierr = vBCvel.end_access(); CHKERRQ(ierr);
    ierr = no_model_mask.end_access(); CHKERRQ(ierr);
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
    modifier = new SIAFD_Regional(grid, *ice, *EC, config);
  } else {
    modifier = new SSBM_Trivial(grid, *ice, *EC, config);
  }

  if (use_ssa_velocity) {
    my_stress_balance = new SSAFD_Regional(grid, *basal, *ice, *EC, config);
  } else {
    my_stress_balance = new SSB_Trivial(grid, *basal, *ice, *EC, config);
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


PetscErrorCode IceRegionalModel::bootstrap_2d(const char *filename) {
  PetscErrorCode ierr;

  ierr = IceModel::bootstrap_2d(filename); CHKERRQ(ierr);

  bool zgwnm;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (!zgwnm) {
    ierr = verbPrintf(2, grid.com, 
      "continuing bootstrapping for IceRegionalModel from file %s\n",filename); CHKERRQ(ierr);
    ierr =  usurfstore.regrid(filename,0.0); CHKERRQ(ierr);
    ierr =    thkstore.regrid(filename,0.0); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceRegionalModel::initFromFile(const char *filename) {
  PetscErrorCode  ierr;

  ierr = IceModel::initFromFile(filename); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing IceRegionalModel from NetCDF file '%s'...\n",
     filename); CHKERRQ(ierr);

  NCTool nc(grid.com, grid.rank);
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  int last_record;  // find index of the last record in the file
  ierr = nc.get_nrecords(last_record); CHKERRQ(ierr);
  last_record -= 1;

  bool nmm_exists;
  ierr = nc.find_variable("no_model_mask", NULL, nmm_exists); CHKERRQ(ierr);
  if (nmm_exists) {
    ierr = verbPrintf(2,grid.com,
	"  reading 'no_model_mask' from %s ...\n",
	filename); CHKERRQ(ierr);
    // note: communication to fill stencil width should occur inside this call
    ierr = no_model_mask.read(filename, last_record); CHKERRQ(ierr);
  }
  
  bool zgwnm, us_exists, ts_exists;
  ierr = PISMOptionsIsSet("-zero_grad_where_no_model", zgwnm); CHKERRQ(ierr);
  if (!zgwnm) {
    ierr = nc.find_variable("usurfstore", NULL, us_exists); CHKERRQ(ierr);
    ierr = nc.find_variable("thkstore", NULL, ts_exists); CHKERRQ(ierr);
    if ((us_exists) && (ts_exists)) {
      ierr = verbPrintf(2,grid.com,
	"  reading 'usurfstore' from %s ...\n",
	filename); CHKERRQ(ierr);
      ierr = usurfstore.read(filename, last_record); CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com,
	"  reading 'thkstore' from %s ...\n",
	filename); CHKERRQ(ierr);
      ierr = thkstore.read(filename, last_record); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com,
	"  PISM ERROR (IceRegionalModel):\n"
	"    'usurfstore' and/or 'thkstore' not present in %s ... ENDING!\n"
	"    (use option -zero_grad_where_no_model)\n",
	filename); CHKERRQ(ierr);
      PISMEnd();
    }    
  } else {
    usurfstore.set(0.0);
    thkstore.set(0.0);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  // at this point may or may not have a filled-in no_model_mask variable;
  //   next try to create it from user option, which overrides input file
  PetscReal stripkm;
  bool nmm_realset;
  ierr = PISMOptionsReal("-no_model_strip", 
                         "width in km of strip near boundary in which modeling is turned off",
			 stripkm, nmm_realset); CHKERRQ(ierr);
  if (nmm_realset) {
    ierr = verbPrintf(2, grid.com,
      "  option -no_model_strip read ...\n"
      "  (re)setting boundary strip width to %.2f km ...\n",
      stripkm); CHKERRQ(ierr);
    ierr = set_no_model_strip(1000.0*stripkm); CHKERRQ(ierr);
  } else {
    if (nmm_exists) { // the o.k. case; we just need to warn if option is given without number
      bool no_model_strip_set;
      ierr = PISMOptionsIsSet("-no_model_strip", no_model_strip_set); CHKERRQ(ierr);
      if (no_model_strip_set) {
        ierr = verbPrintf(2, grid.com,
          "\nPISMO WARNING: option '-no_model_strip' seen with no real value.  Option ignored\n"
          "  because no_model_mask variable was read from input file.  Proceeding ...\n\n");
          CHKERRQ(ierr);
      }
    } else { // bad case: we still don't have a no_model_mask and we have to fail
      ierr = verbPrintf(1, grid.com,
        "\nPISMO ERROR: option '-no_model_strip X' not seen.  No no_model_mask variable\n"
        "  found in input file.  ENDING ...\n\n");
        CHKERRQ(ierr);
      PISMEnd();
    }
  }

  return 0;
}


PetscErrorCode IceRegionalModel::set_vars_from_options() {
  PetscErrorCode ierr;
  bool nmstripSet;
  PetscReal stripkm;

  // base class reads the -boot_file option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "IceRegionalModel", ""); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
                    "* Initializing IceRegionalModel variables ...\n"); CHKERRQ(ierr);
  ierr = PISMOptionsReal("-no_model_strip", 
                         "width in km of strip near boundary in which modeling is turned off",
			 stripkm, nmstripSet);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (nmstripSet) {
    ierr = verbPrintf(2, grid.com,
      "    option -no_model_strip read ... setting boundary strip width to %.2f km\n",
      stripkm); CHKERRQ(ierr);
    ierr = set_no_model_strip(convert(stripkm, "km", "m")); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com,
      "\nPISMO ERROR: option '-no_model_strip X' is REQUIRED if '-i' is not used.\n"
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
  // IceRegionalModel::cell_interface_diffusive_flux() below.
  ierr = no_model_mask.begin_access(); CHKERRQ(ierr);

  ierr = IceModel::massContExplicitStep(); CHKERRQ(ierr);

  ierr = no_model_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes diffusive (SIA) fluxes through interfaces of a computational cell.
/*!
 * This disables diffusive (SIA) flow in the no_model_strip.
 */
PetscErrorCode IceRegionalModel::cell_interface_diffusive_flux(IceModelVec2Stag &Qstag, int i, int j,
                                                               planeStar<PetscScalar> &Q_output) {
  PetscErrorCode ierr;
  if (no_model_mask(i, j) > 0.5) {
    Q_output.e = Q_output.w = Q_output.n = Q_output.s = 0;
  } else {
    ierr = IceModel::cell_interface_diffusive_flux(Qstag, i, j, Q_output); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Computes advective velocities through cell interfaces.
/*!
 * This disables advective (SSA) flow in the no_model_strip.
 */
PetscErrorCode IceRegionalModel::cell_interface_velocities(bool do_part_grid,
                                                           int i, int j,
                                                           planeStar<PetscScalar> &v) {
  PetscErrorCode  ierr;
  if (no_model_mask(i, j) > 0.5) {
    v.e = v.w = v.n = v.s = 0;
  } else {
    ierr = IceModel::cell_interface_velocities(do_part_grid, i, j, v); CHKERRQ(ierr);
  }

  return 0;
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
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

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

