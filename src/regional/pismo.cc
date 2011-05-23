// Copyright (C) 2010, 2011 Ed Bueler, Daniella DellaGiustina and Constantine Khroulev
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
#include "grid.hh"
#include "iceModel.hh"

#include "PCFactory.hh"
#include "PISMAtmosphere.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMStressBalance.hh"
#include "SIAFD.hh"
#include "SSAFD.hh"

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
};

PetscErrorCode SIAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SIAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SIA solver...\n"); CHKERRQ(ierr);

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(1, "no_model_mask is not available");

  return 0;
}

PetscErrorCode SIAFD_Regional::compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y) {
  PetscErrorCode ierr;

  ierr = SIAFD::compute_surface_gradient(h_x, h_y); CHKERRQ(ierr);

  ierr = h_x.begin_access(); CHKERRQ(ierr);
  ierr = h_y.begin_access(); CHKERRQ(ierr);
  ierr = no_model_mask->begin_access(); CHKERRQ(ierr);
  PetscInt GHOSTS = 1;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      if ( ((*no_model_mask)(i,j) > 0.5) && ((*no_model_mask)(i+1,j) > 0.5) ) {
        h_x(i,j,0) = 0.0;
        h_y(i,j,0) = 0.0;
      }
      if ( ((*no_model_mask)(i,j) > 0.5) && ((*no_model_mask)(i,j+1) > 0.5) ) {
        h_x(i,j,1) = 0.0;
        h_y(i,j,1) = 0.0;
      }
    }
  }
  ierr = no_model_mask->end_access(); CHKERRQ(ierr);
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
};

PetscErrorCode SSAFD_Regional::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = SSAFD::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"  using the regional version of the SSA solver...\n"); CHKERRQ(ierr);

  no_model_mask = dynamic_cast<IceModelVec2Int*>(vars.get("no_model_mask"));
  if (no_model_mask == NULL) SETERRQ(1, "no_model_mask is not available");
  
  return 0;
}

PetscErrorCode SSAFD_Regional::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  ierr = SSAFD::compute_driving_stress(result); CHKERRQ(ierr);

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = no_model_mask->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((*no_model_mask)(i,j) > 0.5) {
        result(i,j).u = 0.0;
        result(i,j).v = 0.0;
      }
    }
  }
  ierr = no_model_mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, NCConfigVariable &c, NCConfigVariable &o)
     : IceModel(g,c,o) {};

protected:
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode initFromFile(const char *filename);
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode init_physics();
private:
  IceModelVec2Int   no_model_mask;    
};

PetscErrorCode IceRegionalModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  // stencil width of 2 needed for surfaceGradientSIA() action
  ierr = no_model_mask.create(grid, "no_model_mask", true, 2); CHKERRQ(ierr);
  ierr = no_model_mask.set_attrs("model_state", // ensures that it gets written at the end of the run
    "mask: compute driving stress and surface gradient normally or replace by zero",
    "", ""); CHKERRQ(ierr); // no units and no standard name

  double NMMASK_NORMAL   = 0.0,
         NMMASK_ZERO_OUT = 1.0;
  vector<double> mask_values(2);
  mask_values[0] = NMMASK_NORMAL;
  mask_values[1] = NMMASK_ZERO_OUT;
  ierr = no_model_mask.set_attr("flag_values", mask_values); CHKERRQ(ierr);
  ierr = no_model_mask.set_attr("flag_meanings",
			"normal zero_out_driving_stress_and_surface_gradient");
			CHKERRQ(ierr);
  no_model_mask.output_data_type = NC_BYTE;
  ierr = no_model_mask.set(NMMASK_NORMAL); CHKERRQ(ierr);
  ierr = variables.add(no_model_mask); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode IceRegionalModel::init_physics() {
  PetscErrorCode ierr;

  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  delete stress_balance; // because we delete the old one, at run time we expect
                         //   a second initialization message, from code below 

  // Re-create the stress balance object:
  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_sia = config.get_flag("do_sia");
  
  // We always have SIA "on", but SSA is "on" only if use_ssa_velocity is set.
  // In that case SIA and SSA velocities are always added up (there is no
  // switch saying "do the hybrid").
  ierr = verbPrintf(2,grid.com,
    "  old stress balance and modifier deleted;\n"
    "    replacing with versions following no_model_mask semantics ...\n");
    CHKERRQ(ierr);

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
  } else {
    ierr = verbPrintf(2,grid.com,
	"IceRegionalModel/pismo ERROR: Option '-i' was used but input file %s does not\n"
	"  contain 'no_model_mask' variable.  Executable pismo has no well-defined\n"
	"  semantics without it!  Ending!!\n",
	filename); CHKERRQ(ierr);
    PISMEnd();
  } 

  // at this point we *do* have a no_model_mask variable; now warn user if they
  //   had a -no_model_strip option set; we are going to ignor that option
  bool no_model_strip_set;
  ierr = PISMOptionsIsSet("-no_model_strip", no_model_strip_set); CHKERRQ(ierr);
  if (no_model_strip_set) {
    ierr = PetscPrintf(grid.com,
      "\nPISMO WARNING: option '-no_model_strip' or '-no_model_strip X' seen.  Value X ignored\n"
      "  because no_model_mask variable was read from input file.  Proceeding ...\n\n"); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::set_vars_from_options() {
  PetscErrorCode ierr;
  PetscReal stripkm;
  PetscTruth  nmstripSet;

  // base class reads the -boot_file option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "IceRegionalModel", ""); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
                    "* Initializing IceRegionalModel variables ...\n"); CHKERRQ(ierr);

  bool no_model_strip_set;
  ierr = PISMOptionsIsSet("-no_model_strip", no_model_strip_set); CHKERRQ(ierr);
  if (!no_model_strip_set) {
    ierr = PetscPrintf(grid.com,
                       "\nPISMO ERROR: option '-no_model_strip X' is REQUIRED if '-i' is not used.\n"
                       "   Executable pismo has no well-defined semantics without it!  Ending!!\n\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = PetscOptionsReal("-no_model_strip",
			  "Specifies the no-modeling boundary strip width, in km",
			  "", 20.0,
			  &stripkm, &nmstripSet); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (nmstripSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com,
                      "    option -no_model_strip read ... setting boundary strip width to %.2f km\n",
                      stripkm); CHKERRQ(ierr);
    double strip = 1000.0*stripkm;

    ierr = no_model_mask.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if ((grid.x[i] <= grid.x[0]+strip)
            || (grid.x[i] >= grid.x[grid.Mx-1]-strip)) {
          no_model_mask(i, j) = 1; CHKERRQ(ierr);
        } else if ((grid.y[j] <= grid.y[0]+strip)
                   || (grid.y[j] >= grid.y[grid.My-1]-strip)) {
          no_model_mask(i, j) = 1; CHKERRQ(ierr);
        } else {
          no_model_mask(i, j) = 0; CHKERRQ(ierr);
        }
      }
    }
    ierr = no_model_mask.end_access(); CHKERRQ(ierr);

    ierr = no_model_mask.beginGhostComm(); CHKERRQ(ierr);
    ierr = no_model_mask.endGhostComm(); CHKERRQ(ierr);
  }

  return 0;
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
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
      "  pismo {-i IN.nc|-boot_file IN.nc} [OTHER PISM & PETSc OPTIONS]\n"
      "where:\n"
      "  -i          IN.nc is input file in NetCDF format: contains PISM-written model state\n"
      "  -boot_file  IN.nc is input file in NetCDF format: contains a few fields, from which\n"
      "              heuristics will build initial model state\n"
      "notes:\n"
      "  * one of -i or -boot_file is required\n"
      "  * if -boot_file is used then also '-Mx A -My B -Mz C -Lz D' are required\n";
    if ((!iset) && (!bfset)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismo", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismo", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    // initialize the ice dynamics model
    IceGrid g(com, rank, size, config);
    IceRegionalModel m(g, config, overrides);

    // initialize boundary models
    PAFactory pa(g, config);
    PISMAtmosphereModel *atmosphere;

    PSFactory ps(g, config);
    PISMSurfaceModel *surface;

    POFactory po(g, config);
    PISMOceanModel *ocean;

    ierr = PetscOptionsBegin(com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);
    pa.create(atmosphere);
    ps.create(surface);
    po.create(ocean);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    surface->attach_atmosphere_model(atmosphere);

    m.attach_ocean_model(ocean);
    m.attach_surface_model(surface);
    ierr = m.setExecName("pismo"); CHKERRQ(ierr);

    ierr = m.init(); CHKERRQ(ierr);

    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);
    // provide a default output file name if no -o option is given.
    ierr = m.writeFiles("unnamed_regional.nc"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

