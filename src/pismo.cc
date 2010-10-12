// Copyright (C) 2010 Ed Bueler and Daniella DellaGiustina
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
#include "base/grid.hh"
#include "base/iceModel.hh"

#include "coupler/PCFactory.hh"
#include "coupler/PISMAtmosphere.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/PISMOcean.hh"

class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid &g, NCConfigVariable &config, NCConfigVariable &overrides)
     : IceModel(g,config,overrides) {};

protected:
  virtual PetscErrorCode set_vars_from_options();
  virtual PetscErrorCode initFromFile(const char *filename);
  virtual PetscErrorCode write_extra_fields(const char* filename);
  virtual PetscErrorCode createVecs();
  virtual PetscErrorCode computeDrivingStress(IceModelVec2S &vtaudx, IceModelVec2S &vtaudy);
  virtual PetscErrorCode surfaceGradientSIA();

private:
    IceModelVec2S   no_model_mask;    
};


PetscErrorCode IceRegionalModel::createVecs() {
  PetscErrorCode ierr;

  ierr = no_model_mask.create(grid, "no_model_mask", false); CHKERRQ(ierr);
  ierr = no_model_mask.set_attrs("internal",
                            "mask specifying whether to model the ice sheet (=0), or hold to boundary values already set by input file [or externally in future] (=1)",
                            "", ""); CHKERRQ(ierr); // no units and no standard name
  ierr = no_model_mask.set(0.0); CHKERRQ(ierr);    // set to no such strip of boundary values

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
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
  ierr = nc.get_dim_length("t", &last_record); CHKERRQ(ierr);
  last_record -= 1;
  bool nmm_exists;
  ierr = nc.find_variable("no_model_mask", NULL, nmm_exists); CHKERRQ(ierr);

  if (nmm_exists) {
    ierr = verbPrintf(2,grid.com,
	"  reading 'no_model_mask' from %s ...\n",
	filename); CHKERRQ(ierr);
    ierr = no_model_mask.read(filename, last_record); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com,
	"IceRegionalModel WARNING: input file %s does not contain 'no_model_mask'\n"
	"  variable.  Keeping its value as identically zero (= everywhere modeled).\n",
	filename); CHKERRQ(ierr);
  } 

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::set_vars_from_options() {
  PetscErrorCode ierr;

  // base class reads the -boot_file option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "IceRegionalModel", ""); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com, 
     "* Initializing IceRegionalModel variables ...\n"); CHKERRQ(ierr);

  PetscReal stripkm;
  PetscTruth  nmstripSet;
  ierr = PetscOptionsReal("-no_model_strip",
			  "Specifies the no-modeling boundary strip width, in km",
			  "", 20.0,
			  &stripkm, &nmstripSet); CHKERRQ(ierr);

  if (nmstripSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com,
       "    option -no_model_strip read ... setting boundary strip width to %.2f km\n",
       stripkm); CHKERRQ(ierr);
    double *x_coords, *y_coords, strip = 1000.0*stripkm;
    ierr = grid.compute_horizontal_coordinates(x_coords, y_coords); CHKERRQ(ierr);

    ierr = no_model_mask.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if ((x_coords[i] <= x_coords[0]+strip)
            || (x_coords[i] >= x_coords[grid.Mx-1]-strip)) {
          no_model_mask(i, j) = 1; CHKERRQ(ierr);
        } else if ((y_coords[j] <= y_coords[0]+strip)
                   || (y_coords[j] >= y_coords[grid.My-1]-strip)) {
          no_model_mask(i, j) = 1; CHKERRQ(ierr);
        } else {
          no_model_mask(i, j) = 0; CHKERRQ(ierr);
        }
      }
    }
    ierr = no_model_mask.end_access(); CHKERRQ(ierr);

    delete[] x_coords;  delete[] y_coords;
  }

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::write_extra_fields(const char* filename) {
  PetscErrorCode ierr;
  ierr = IceModel::write_extra_fields(filename); CHKERRQ(ierr);

  ierr = no_model_mask.write(filename); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceRegionalModel::computeDrivingStress(
                    IceModelVec2S &vtaudx, IceModelVec2S &vtaudy) {
  PetscErrorCode ierr;

  ierr = IceModel::computeDrivingStress(vtaudx, vtaudy); CHKERRQ(ierr);
                    
  ierr = vtaudx.begin_access(); CHKERRQ(ierr);
  ierr = vtaudy.begin_access(); CHKERRQ(ierr);
  ierr = no_model_mask.begin_access();  CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (no_model_mask(i,j) > 0.5) {
        vtaudx(i,j) = 0.0;
        vtaudy(i,j) = 0.0;
      }
    }
  }
  ierr = no_model_mask.end_access(); CHKERRQ(ierr);
  ierr = vtaudx.end_access(); CHKERRQ(ierr);
  ierr = vtaudy.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceRegionalModel::surfaceGradientSIA() {
  PetscErrorCode  ierr;

  ierr = IceModel::surfaceGradientSIA(); CHKERRQ(ierr);

  PetscScalar **h_x[2], **h_y[2];
  ierr = vWork2d[0].get_array(h_x[0]); CHKERRQ(ierr);
  ierr = vWork2d[1].get_array(h_x[1]); CHKERRQ(ierr);
  ierr = vWork2d[2].get_array(h_y[0]); CHKERRQ(ierr);
  ierr = vWork2d[3].get_array(h_y[1]); CHKERRQ(ierr);
  ierr = no_model_mask.begin_access();  CHKERRQ(ierr);

  // FIXME:  because no_model_mask does not have stencil width, we are not doing
  //   ghosts here, as is done in iMsia.cc
  for (PetscInt o=0; o<2; o++) {
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        if (no_model_mask(i,j) > 0.5) {
          h_x[o][i][j] = 0.0;
          h_y[o][i][j] = 0.0;
        }
      }
    }
  }

  ierr = no_model_mask.end_access();  CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[2].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[3].end_access(); CHKERRQ(ierr);

  // FIXME:  probably this communication is necessary for now; probably no_model_mask should
  //   have ghosts of stencil width 1 so that the communication can be removed
  //   and the loop above can do local update of ghosts like the loop in iMsia.cc
  for (PetscInt k=0; k<4; k++) {
    ierr = vWork2d[k].beginGhostComm(); CHKERRQ(ierr);
    ierr = vWork2d[k].endGhostComm(); CHKERRQ(ierr);
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

    ierr = check_old_option_and_stop(com, "-boot_from", "-boot_file"); CHKERRQ(ierr); 

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
    if ((iset == PETSC_FALSE) && (bfset == PETSC_FALSE)) {
      ierr = PetscPrintf(com,
         "\nPISM ERROR: one of options -i,-boot_file is required\n\n"); CHKERRQ(ierr);
      ierr = show_usage_and_quit(com, "pismo", usage.c_str()); CHKERRQ(ierr);
    } else {
      vector<string> required;  required.clear();
      ierr = show_usage_check_req_opts(com, "pismo", required, usage.c_str()); CHKERRQ(ierr);
    }

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid g(com, rank, size, config);
//    IceModel m(g, config, overrides);
    IceRegionalModel m(g, config, overrides);

    // Initialize boundary models:
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
