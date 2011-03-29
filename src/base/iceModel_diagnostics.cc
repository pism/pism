// Copyright (C) 2010, 2011 Constantine Khroulev
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

#include "iceModel_diagnostics.hh"
#include "PISMDiagnostic.hh"


PetscErrorCode IceModel::init_diagnostics() {

  // Add IceModel diagnostics:
  diagnostics["bwp"]              = new IceModel_bwp(this, grid, variables);
  diagnostics["cts"]              = new IceModel_cts(this, grid, variables);
  diagnostics["dhdt"]             = new IceModel_dhdt(this, grid, variables);
  diagnostics["enthalpybase"]     = new IceModel_enthalpybase(this, grid, variables);
  diagnostics["enthalpysurf"]     = new IceModel_enthalpysurf(this, grid, variables);
  diagnostics["hardav"]           = new IceModel_hardav(this, grid, variables);
  diagnostics["liqfrac"]          = new IceModel_liqfrac(this, grid, variables);
  diagnostics["proc_ice_area"]    = new IceModel_proc_ice_area(this, grid, variables);
  diagnostics["rank"]             = new IceModel_rank(this, grid, variables);
  diagnostics["temp"]             = new IceModel_temp(this, grid, variables);
  diagnostics["temp_pa"]          = new IceModel_temp_pa(this, grid, variables);
  diagnostics["tempbase"]         = new IceModel_tempbase(this, grid, variables);
  diagnostics["tempicethk"]       = new IceModel_tempicethk(this, grid, variables);
  diagnostics["tempicethk_basal"] = new IceModel_tempicethk_basal(this, grid, variables);
  diagnostics["temppabase"]       = new IceModel_temppabase(this, grid, variables);
  diagnostics["tempsurf"]         = new IceModel_tempsurf(this, grid, variables);
  diagnostics["new_mask"]         = new IceModel_new_mask(this, grid, variables);

  // Get diagnostics supported by the stress balance object:
  stress_balance->get_diagnostics(diagnostics);

  // Get diagnostics supported by the surface model:
  surface->get_diagnostics(diagnostics);

  // Get diagnostics supported by the ocean model:
  ocean->get_diagnostics(diagnostics);

  // Get diagnostics supported by the bed deformation model:
  if (beddef) {
    beddef->get_diagnostics(diagnostics);
  }

  int threshold = 5;
  if (getVerbosityLevel() >= threshold) {
    verbPrintf(threshold, grid.com, " *** Available diagnostic quantities:\n");

    map<string, PISMDiagnostic*>::iterator j = diagnostics.begin();
    while (j != diagnostics.end()) {
      string name = j->first;
      PISMDiagnostic *diag = j->second;

      int N = diag->get_nvars();
      verbPrintf(threshold, grid.com, " ** %s\n", name.c_str());

      for (int k = 0; k < N; ++k) {
        NCSpatialVariable *var = diag->get_metadata(k);

        string long_name = var->get_string("long_name");

        verbPrintf(threshold, grid.com, " * %s\n", long_name.c_str());
      }
      ++j;
    }
  }

  return 0;
}

IceModel_hardav::IceModel_hardav(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("hardav", grid);

  const PetscScalar power = 1.0 / model->ice->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  vars[0].set("valid_min", 0);
  vars[0].set("_FillValue", -0.01);
}

//! \brief Computes vertically-averaged ice hardness.
PetscErrorCode IceModel_hardav::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  const PetscScalar fillval = -0.01;
  PetscScalar *Eij; // columns of enthalpy values

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "hardav", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = model->Enth3.getInternalColumn(i,j,&Eij); CHKERRQ(ierr);
      const PetscScalar H = model->vH(i,j);
      if (H > 0.0) {
        (*result)(i,j) = model->ice->averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                                grid.zlevels.data(), Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


IceModel_rank::IceModel_rank(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("rank", grid);
  
  set_attrs("processor rank", "", "", "", 0);
  vars[0].time_independent = true;
}

PetscErrorCode IceModel_rank::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "rank", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      (*result)(i,j) = grid.rank;
  ierr = result->end_access();

  output = result;
  return 0;
}

//! Computes the subglacial (basal) water pressure
IceModel_bwp::IceModel_bwp(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("bwp", grid);
  set_attrs("subglacial (pore) water pressure", "", "Pa", "Pa", 0);
  vars[0].set("_FillValue", -0.01);
  vars[0].set("valid_min", 0);
}

/*!
  \f[p_w = \alpha\, \frac{w}{w_{\text{max}}}\, \rho\, g\, H,\f]
  where 

  - \f$\alpha\f$ is the till pore water fraction (till_pw_fraction),
  - \f$w\f$ is the effective thickness of subglacial melt water (bwat)
  - \f$w_{\text{max}}\f$ is the maximum allowed value for \f$w\f$ (hmelt_max),
  - \f$\rho\f$ is the ice density (ice_density)
  - \f$H\f$ is the ice thickness (thk)

Result is set to invalid (_FillValue) where the ice is floating, there being
no meaning to the above calculation.
 */
PetscErrorCode IceModel_bwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  const PetscScalar
    alpha     = model->config.get("till_pw_fraction"),
    wmax      = model->config.get("hmelt_max"),
    fillval   = -0.01;

  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = model->vHmelt.begin_access(); CHKERRQ(ierr);
  ierr = model->vbmr.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (model->vH(i,j) > 0.0) {
        (*result)(i,j) = model->getBasalWaterPressure(model->vH(i,j), // FIXME task #7297
                                                   model->vHmelt(i,j),
                                                   model->vbmr(i,j),
                                                   alpha, wmax);
      } else { // put negative value below valid range
        (*result)(i,j) = fillval;
      }
    }
  }
  ierr = model->vH.end_access(); CHKERRQ(ierr);
  ierr = model->vHmelt.end_access(); CHKERRQ(ierr);
  ierr = model->vbmr.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = model->vMask.fill_where_floating(*result, fillval); CHKERRQ(ierr);
  
  output = result;
  return 0;
}


IceModel_cts::IceModel_cts(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_3d("cts", grid, g.zlevels);
  
  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

PetscErrorCode IceModel_cts::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "cts", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->setCTSFromEnthalpy(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

//! \brief Computes the rate of change of ice surface elevation.
IceModel_dhdt::IceModel_dhdt(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("dhdt", grid);
  
  set_attrs("rate of change of surface elevation", "",
            "m s-1", "m year-1", 0);
}

PetscErrorCode IceModel_dhdt::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "dhdt", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->vdHdt.begin_access(); CHKERRQ(ierr);
  ierr = model->vuplift.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      (*result)(i,j) = model->vdHdt(i,j) + model->vuplift(i,j);
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vuplift.end_access(); CHKERRQ(ierr);
  ierr = model->vdHdt.end_access(); CHKERRQ(ierr);
  
  output = result;
  return 0;
}

IceModel_proc_ice_area::IceModel_proc_ice_area(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("proc_ice_area", grid);
  
  set_attrs("number of cells containing ice in a processor's domain", "",
            "", "", 0);
  vars[0].time_independent = true;
}

PetscErrorCode IceModel_proc_ice_area::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *thickness;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "proc_ice_area", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PetscInt ice_filled_cells = 0;

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      if ((*thickness)(i,j) > 0) {
	ice_filled_cells += 1;
      }
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i)
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j)
      (*result)(i,j) = ice_filled_cells;
  ierr = result->end_access();
  
  output = result;
  return 0;
}


IceModel_temp::IceModel_temp(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_3d("temp", grid, g.zlevels);
  
  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  vars[0].set("valid_min", 0);
}

PetscErrorCode IceModel_temp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);
  
  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "temp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");

  PetscScalar *Tij, *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result->getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = (*thickness)(i,j) - grid.zlevels[k];
        ierr = model->EC->getAbsTemp(Enthij[k],
                                     model->EC->getPressureFromDepth(depth),
                                     Tij[k]); 
        if (ierr) {
          PetscPrintf(grid.com,
            "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d,k=%d\n\n",
            i,j,k);
        }
        CHKERRQ(ierr);
      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);
  
  output = result;
  return 0;
}


IceModel_temp_pa::IceModel_temp_pa(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_3d("temp_pa", grid, g.zlevels);
  
  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
}

PetscErrorCode IceModel_temp_pa::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  bool cold_mode = model->config.get_flag("do_cold_ice_methods");
  PetscReal triple_point_temp = model->config.get("water_triple_point_temperature");

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "temp_pa", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");

  PetscScalar *Tij, *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result->getInternalColumn(i,j,&Tij); CHKERRQ(ierr);
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);
      for (PetscInt k=0; k<grid.Mz; ++k) {
        const PetscScalar depth = (*thickness)(i,j) - grid.zlevels[k],
          p = model->EC->getPressureFromDepth(depth);
        ierr = model->EC->getPATemp(Enthij[k], p, Tij[k]); 
        if (ierr) {
          PetscPrintf(grid.com,
            "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d,k=%d\n\n",
            i,j,k);
        }
        CHKERRQ(ierr);

        if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
                         // is 273.15
          if ( model->EC->isTemperate(Enthij[k],p) && ((*thickness)(i,j) > 0)) {
            Tij[k] = triple_point_temp;
          }
        }

      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->shift(-triple_point_temp); CHKERRQ(ierr); 

  output = result;
  return 0;
}

IceModel_temppabase::IceModel_temppabase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("temppabase", grid);
  
  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "degrees Celsius", "degrees Celsius", 0);
}

PetscErrorCode IceModel_temppabase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  bool cold_mode = model->config.get_flag("do_cold_ice_methods");
  PetscReal triple_point_temp = model->config.get("water_triple_point_temperature");
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "temp_pa_base", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  IceModelVec3 *enthalpy;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(variables.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");

  PetscScalar *Enthij; // columns of these values
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&Enthij); CHKERRQ(ierr);

      const PetscScalar depth = (*thickness)(i,j),
        p = model->EC->getPressureFromDepth(depth);
      ierr = model->EC->getPATemp(Enthij[0], p,
                                  (*result)(i,j)); 
      if (ierr) {
        PetscPrintf(grid.com,
                    "\n\nEnthalpyConverter.getAbsTemp() error at i=%d,j=%d\n\n",
                    i,j);
      }
      CHKERRQ(ierr);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
                       // is 273.15
        if ( model->EC->isTemperate(Enthij[0],p) && ((*thickness)(i,j) > 0)) {
          (*result)(i,j) = triple_point_temp;
        }
      }
    }
  }
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->shift(-triple_point_temp); CHKERRQ(ierr); 
  
  output = result;
  return 0;
}

IceModel_enthalpysurf::IceModel_enthalpysurf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("enthalpysurf", grid);
  
  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set("_FillValue", -0.01);
}

PetscErrorCode IceModel_enthalpysurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "enthalpysurf", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;

  // compute levels corresponding to 1 m below the ice surface:

  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      (*result)(i,j) = PetscMax(model->vH(i,j) - 1.0, 0.0);
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = model->Enth3.getSurfaceValues(*result, *result); CHKERRQ(ierr);  // z=0 slice

  ierr = result->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (model->vH(i,j) <= 1.0)
	(*result)(i,j) = fill_value;
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);

  
  output = result;
  return 0;
}

IceModel_enthalpybase::IceModel_enthalpybase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("enthalpybase", grid);
  
  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set("_FillValue", -0.01);
}

PetscErrorCode IceModel_enthalpybase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "enthalpybase", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->Enth3.getHorSlice(*result, 0.0); CHKERRQ(ierr);  // z=0 slice

  ierr = result->mask_by(model->vH, -0.01); CHKERRQ(ierr); 

  output = result;
  return 0;
}


IceModel_tempbase::IceModel_tempbase(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("tempbase", grid);
  
  set_attrs("ice temperature at the base of ice", "",
            "K", "K", 0);
  vars[0].set("_FillValue", -0.01);
}

PetscErrorCode IceModel_tempbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result, *thickness;

  IceModel_enthalpybase enth(model, grid, variables);

  ierr = enth.compute(output); CHKERRQ(ierr);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) SETERRQ(1, "dynamic_cast failure");

  // result contains basal enthalpy; note that it is allocated by
  // enth.compute().

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      PetscReal depth = (*thickness)(i,j),
        pressure = model->EC->getPressureFromDepth(depth);
      if (depth > 0) {
        ierr = model->EC->getAbsTemp((*result)(i,j),
                                     pressure,
                                     (*result)(i,j));
        CHKERRQ(ierr);
      } else {
        (*result)(i,j) = -0.01;
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  output = result;
  return 0;
}

IceModel_tempsurf::IceModel_tempsurf(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("tempsurf", grid);
  
  set_attrs("ice temperature at 1m below the ice surface", "",
            "K", "K", 0);
  vars[0].set("_FillValue", -0.01);
}

PetscErrorCode IceModel_tempsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result, *thickness;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  IceModel_enthalpysurf enth(model, grid, variables);

  ierr = enth.compute(output); CHKERRQ(ierr);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) SETERRQ(1, "dynamic_cast failure");

  // result contains surface enthalpy; note that it is allocated by
  // enth.compute().

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  PetscReal depth = 1.0,
    pressure = model->EC->getPressureFromDepth(depth);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*thickness)(i,j) > 1) {
        ierr = model->EC->getAbsTemp((*result)(i,j),
                                     pressure,
                                     (*result)(i,j)); CHKERRQ(ierr);
      } else {
        (*result)(i,j) = -0.01;
      }
    }
  }
  
  ierr = thickness->end_access(); CHKERRQ(ierr);  
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  output = result;
  return 0;
}


IceModel_liqfrac::IceModel_liqfrac(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_3d("liqfrac", grid, g.zlevels);
  
  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
}

PetscErrorCode IceModel_liqfrac::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.zlevels);

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "liqfrac", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  bool cold_mode = model->config.get_flag("do_cold_ice_methods");

  if (cold_mode) {
    ierr = result->set(0.0); CHKERRQ(ierr);
  } else {
    ierr = model->compute_liquid_water_fraction(model->Enth3, *result); CHKERRQ(ierr);
  }

  output = result;
  return 0;
}

IceModel_tempicethk::IceModel_tempicethk(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  PetscScalar fill_value = -0.01;

  // set metadata:
  vars[0].init_2d("tempicethk", grid);
  
  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  vars[0].set("_FillValue", fill_value);
}

PetscErrorCode IceModel_tempicethk::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "tempicethk", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PetscScalar *Enth;

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (model->vH(i,j) > 0.) {
	ierr = model->Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
	PetscScalar tithk = 0.;
	const PetscInt ks = grid.kBelowHeight(model->vH(i,j));
        
	for (PetscInt k=0; k<ks; ++k) { // FIXME task #7297
          PetscReal pressure = model->EC->getPressureFromDepth(model->vH(i,j) - grid.zlevels[k]);

	  if (model->EC->isTemperate(Enth[k], pressure)) {
	    tithk += grid.zlevels[k+1] - grid.zlevels[k];
	  }
	}

        PetscReal pressure = model->EC->getPressureFromDepth(model->vH(i,j) - grid.zlevels[ks]);
	if (model->EC->isTemperate(Enth[ks], pressure)) {
	  tithk += model->vH(i,j) - grid.zlevels[ks];
	}

	(*result)(i,j) = tithk;
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  PetscScalar fill_value = -0.01;
  ierr = result->mask_by(model->vH, fill_value); CHKERRQ(ierr);

  output = result;
  return 0;
}

IceModel_tempicethk_basal::IceModel_tempicethk_basal(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  PetscScalar fill_value = -0.01;
  
  // set metadata:
  vars[0].init_2d("tempicethk_basal", grid);
  
  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  vars[0].set("_FillValue", fill_value);
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
PetscErrorCode IceModel_tempicethk_basal::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "tempicethk_basal", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PetscScalar *Enth, fill_value = -0.01;
  EnthalpyConverter *EC = model->EC;

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = model->Enth3.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      PetscReal thk = model->vH(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (thk < 0.1)
        continue;

      ierr = model->Enth3.getInternalColumn(i,j,&Enth); CHKERRQ(ierr);
      PetscReal pressure;
      PetscInt ks = grid.kBelowHeight(thk),
        k = 0;

      while (k <= ks) {         // FIXME task #7297
        pressure = EC->getPressureFromDepth(thk - grid.zlevels[k]);

        if (EC->isTemperate(Enth[k],pressure))
          k++;
        else
          break;
      }
      // after this loop 'pressure' is equal to the pressure at the first level
      // that is cold

      // no temperate ice at all; go to the next grid point
      if (k == 0) {
        (*result)(i,j) = 0;
        continue;
      }
      
      // the whole column is temperate (except, possibly, some ice between
      // zlevels[ks] and the total thickness; we ignore it)
      if (k == ks + 1) {
        (*result)(i,j) = grid.zlevels[ks];
        continue;
      }

      PetscReal 
        pressure_0 = EC->getPressureFromDepth(thk - grid.zlevels[k-1]),
        dz = grid.zlevels[k] - grid.zlevels[k-1],
        slope1 = (Enth[k] - Enth[k-1]) / dz,
        slope2 = (EC->getEnthalpyCTS(pressure) - EC->getEnthalpyCTS(pressure_0)) / dz;
      
      if (slope1 != slope2) {
        (*result)(i,j) = grid.zlevels[k-1] +
          (EC->getEnthalpyCTS(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i,j) = PetscMax((*result)(i,j), grid.zlevels[k-1]);
        (*result)(i,j) = PetscMin((*result)(i,j), grid.zlevels[k]);
      } else {
        SETERRQ4(1, "This should never happen: (i=%d, j=%d, k=%d, ks=%d)\n",
                    i, j, k, ks);
      }
    }
  }
  ierr = model->Enth3.end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);

  ierr = result->mask_by(model->vH, fill_value); CHKERRQ(ierr);
  
  output = result;
  return 0;
}


IceModel_new_mask::IceModel_new_mask(IceModel *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<IceModel>(m, g, my_vars) {
  
  // set metadata:
  vars[0].init_2d("new_mask", grid);
  
  vector<double> values(3);
  values[0] = 1;
  values[1] = 2;
  values[2] = 4;

  vars[0].doubles["flag_masks"] = values;
  vars[0].set_string("flag_meanings", "ice_filled_cell interior_cell floating_or_ocean_cell");
  vars[0].set_string("long_name", "new full/ice-free, grounded/ocean, interior/boundary binary mask");
}

PetscErrorCode IceModel_new_mask::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  
  IceModelVec2Mask *result = new IceModelVec2Mask;
  ierr = result->create(grid, "new_mask", true); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  if (model->ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal currentSeaLevel;
  ierr = model->ocean->sea_level_elevation(grid.year, model->dt / secpera, currentSeaLevel); CHKERRQ(ierr);

  double ocean_rho = model->config.get("sea_water_density");

  ierr = result->begin_access(); CHKERRQ(ierr);

  // clear the ghosts
  ierr = result->set(0.0); CHKERRQ(ierr); 

  ierr = model->vH.begin_access(); CHKERRQ(ierr);
  ierr = model->vbed.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      const PetscScalar hgrounded = model->vbed(i,j) + model->vH(i,j), // FIXME task #7297
        hfloating = currentSeaLevel + (1.0 - model->ice->rho/ocean_rho) * model->vH(i,j);

      bool is_floating = hfloating > hgrounded + 1.0,
        // note: the following implies that ice-free cells with bed evelation
        // exactly at sea level are considered grounded
        has_ice = model->vH(i,j) > 0.01;

      int mask_value = 0;
      if (is_floating)
        mask_value = mask_value | FLAG_IS_OCEAN;

      if (has_ice)
        mask_value = mask_value | FLAG_IS_FULL;

      (*result)(i,j) = mask_value;
    }
  }
  ierr = model->vbed.end_access(); CHKERRQ(ierr);
  ierr = model->vH.end_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {

      int mask_ij = result->value(i,j),
        mask_e = result->value(i + 1, j),
        mask_w = result->value(i - 1, j),
        mask_n = result->value(i, j + 1),
        mask_s = result->value(i, j - 1);

      // check if a point is in the interior of a region covered with ice:
      bool ice_filled_interior =
        (mask_ij & FLAG_IS_FULL) &&
        (mask_e  & FLAG_IS_FULL) &&
        (mask_w  & FLAG_IS_FULL) &&
        (mask_s  & FLAG_IS_FULL) &&
        (mask_n  & FLAG_IS_FULL);

      // check if a point is in the interior of an ice-free region:
      bool ice_free_interior =
        (! (mask_ij & FLAG_IS_FULL) ) &&
        (! (mask_e  & FLAG_IS_FULL) ) &&
        (! (mask_w  & FLAG_IS_FULL) ) &&
        (! (mask_s  & FLAG_IS_FULL) ) &&
        (! (mask_n  & FLAG_IS_FULL) );
      
      if (ice_filled_interior || ice_free_interior)
        mask_ij = mask_ij | FLAG_IS_INTERIOR;

      (*result)(i,j) = mask_ij;
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}
