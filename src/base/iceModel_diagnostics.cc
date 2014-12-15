// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "pism_options.hh"
#include "iceModel_diagnostics.hh"
#include "PISMDiagnostic.hh"
#include "Mask.hh"

#include "PISMBedDef.hh"
#include "PISMYieldStress.hh"
#include "PISMHydrology.hh"
#include "PISMStressBalance.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "enthalpyConverter.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "bedrockThermalUnit.hh"
#include "error_handling.hh"
#include "iceModelVec3Custom.hh"

namespace pism {

void IceModel::init_diagnostics() {

  // Add IceModel diagnostics:
  diagnostics["cts"]              = new IceModel_cts(this, grid, variables);
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
  diagnostics["dHdt"]             = new IceModel_dHdt(this, grid, variables);

  if (flux_divergence.was_created()) {
    diagnostics["flux_divergence"] = new IceModel_flux_divergence(this, grid, variables);
  }

  if (climatic_mass_balance_cumulative.was_created()) {
    diagnostics["climatic_mass_balance_cumulative"] = new IceModel_climatic_mass_balance_cumulative(this, grid, variables);
  }

  if (nonneg_flux_2D_cumulative.was_created()) {
    diagnostics["nonneg_flux_cumulative"] = new IceModel_nonneg_flux_2D_cumulative(this, grid, variables);
  }

  if (grounded_basal_flux_2D_cumulative.was_created()) {
    diagnostics["grounded_basal_flux_cumulative"] = new IceModel_grounded_basal_flux_2D_cumulative(this, grid, variables);
  }

  if (floating_basal_flux_2D_cumulative.was_created()) {
    diagnostics["floating_basal_flux_cumulative"] = new IceModel_floating_basal_flux_2D_cumulative(this, grid, variables);
  }

  if (discharge_flux_2D_cumulative.was_created()) {
    diagnostics["discharge_flux_cumulative"] = new IceModel_discharge_flux_2D_cumulative(this, grid, variables);
  }

#if (PISM_USE_PROJ4==1)
  if (global_attributes.has_attribute("proj4")) {
    std::string proj4 = global_attributes.get_string("proj4");
    diagnostics["lat_bnds"] = new IceModel_lat_lon_bounds(this, grid, variables, "lat", proj4);
    diagnostics["lon_bnds"] = new IceModel_lat_lon_bounds(this, grid, variables, "lon", proj4);
  }
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

  ts_diagnostics["ivol"]          = new IceModel_ivol(this, grid, variables);
  ts_diagnostics["slvol"]         = new IceModel_slvol(this, grid, variables);
  ts_diagnostics["divoldt"]       = new IceModel_divoldt(this, grid, variables);
  ts_diagnostics["iarea"]         = new IceModel_iarea(this, grid, variables);
  ts_diagnostics["imass"]         = new IceModel_imass(this, grid, variables);
  ts_diagnostics["dimassdt"]      = new IceModel_dimassdt(this, grid, variables);
  ts_diagnostics["ivoltemp"]      = new IceModel_ivoltemp(this, grid, variables);
  ts_diagnostics["ivolcold"]      = new IceModel_ivolcold(this, grid, variables);
  ts_diagnostics["ivolg"]         = new IceModel_ivolg(this, grid, variables);
  ts_diagnostics["ivolf"]         = new IceModel_ivolf(this, grid, variables);
  ts_diagnostics["iareatemp"]     = new IceModel_iareatemp(this, grid, variables);
  ts_diagnostics["iareacold"]     = new IceModel_iareacold(this, grid, variables);
  ts_diagnostics["iareag"]        = new IceModel_iareag(this, grid, variables);
  ts_diagnostics["iareaf"]        = new IceModel_iareaf(this, grid, variables);
  ts_diagnostics["dt"]            = new IceModel_dt(this, grid, variables);
  ts_diagnostics["max_diffusivity"] = new IceModel_max_diffusivity(this, grid, variables);
  ts_diagnostics["ienthalpy"]     = new IceModel_ienthalpy(this, grid, variables);
  ts_diagnostics["max_hor_vel"]   = new IceModel_max_hor_vel(this, grid, variables);

  ts_diagnostics["surface_ice_flux"]   = new IceModel_surface_flux(this, grid, variables);
  ts_diagnostics["surface_ice_flux_cumulative"]   = new IceModel_surface_flux_cumulative(this, grid, variables);
  ts_diagnostics["grounded_basal_ice_flux"]     = new IceModel_grounded_basal_flux(this, grid, variables);
  ts_diagnostics["grounded_basal_ice_flux_cumulative"]     = new IceModel_grounded_basal_flux_cumulative(this, grid, variables);
  ts_diagnostics["sub_shelf_ice_flux"] = new IceModel_sub_shelf_flux(this, grid, variables);
  ts_diagnostics["sub_shelf_ice_flux_cumulative"] = new IceModel_sub_shelf_flux_cumulative(this, grid, variables);
  ts_diagnostics["nonneg_rule_flux"]   = new IceModel_nonneg_flux(this, grid, variables);
  ts_diagnostics["nonneg_rule_flux_cumulative"]   = new IceModel_nonneg_flux_cumulative(this, grid, variables);
  ts_diagnostics["discharge_flux"]    = new IceModel_discharge_flux(this, grid, variables);
  ts_diagnostics["discharge_flux_cumulative"]    = new IceModel_discharge_flux_cumulative(this, grid, variables);
  ts_diagnostics["H_to_Href_flux"] = new IceModel_H_to_Href_flux(this, grid, variables);
  ts_diagnostics["Href_to_H_flux"] = new IceModel_Href_to_H_flux(this, grid, variables);
  ts_diagnostics["sum_divQ_flux"]  = new IceModel_sum_divQ_flux(this, grid, variables);

  // Get diagnostics supported by the stress balance object:
  if (stress_balance != NULL) {
    stress_balance->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the surface model:
  if (surface != NULL) {
    surface->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the ocean model:
  if (ocean != NULL) {
    ocean->get_diagnostics(diagnostics, ts_diagnostics);
  }

  // Get diagnostics supported by the bed deformation model:
  if (beddef != NULL) {
    beddef->get_diagnostics(diagnostics, ts_diagnostics);
  }

  if (basal_yield_stress_model != NULL) {
    basal_yield_stress_model->get_diagnostics(diagnostics, ts_diagnostics);
  }

  if (subglacial_hydrology != NULL) {
    subglacial_hydrology->get_diagnostics(diagnostics, ts_diagnostics);
  }
}

void IceModel::list_diagnostics() {

  PetscPrintf(grid.com, "\n");

  // quantities with dedicated storage
  {
    std::set<std::string> list = variables.keys();

    if (beddef != NULL) {
      beddef->add_vars_to_output("big", list);
    }

    if (btu != NULL) {
      btu->add_vars_to_output("big", list);
    }

    if (basal_yield_stress_model != NULL) {
      basal_yield_stress_model->add_vars_to_output("big", list);
    }

    if (subglacial_hydrology != NULL) {
      subglacial_hydrology->add_vars_to_output("big", list);
    }

    if (stress_balance != NULL) {
      stress_balance->add_vars_to_output("big", list);
    }

    if (ocean != NULL) {
      ocean->add_vars_to_output("big", list);
    }

    if (surface != NULL) {
      surface->add_vars_to_output("big", list);
    }

    for (unsigned int d = 3; d > 1; --d) {

      PetscPrintf(grid.com,
                  "======== Available %dD quantities with dedicated storage ========\n",
                  d);

      std::set<std::string>::iterator j = list.begin();
      while(j != list.end()) {
        IceModelVec *v = variables.get(*j);

        if (v != NULL && v->get_ndims() == d) {
          NCSpatialVariable var = v->metadata();

          std::string
            name                = var.get_name(),
            units               = var.get_string("units"),
            glaciological_units = var.get_string("glaciological_units"),
            long_name           = var.get_string("long_name");

          if (glaciological_units.empty() == false) {
            units = glaciological_units;
          }

            PetscPrintf(grid.com,
                        "   Name: %s [%s]\n"
                        "       - %s\n\n", name.c_str(), units.c_str(), long_name.c_str());
        }

        ++j;
      }
    }

  } // done with quantities with dedicated storage

  // 2D and 3D diagnostics
  for (unsigned int d = 3; d > 1; --d) {

    PetscPrintf(grid.com,
                "======== Available %dD diagnostic quantities ========\n",
                d);

    std::map<std::string, Diagnostic*>::iterator j = diagnostics.begin();
    while (j != diagnostics.end()) {
      Diagnostic *diag = j->second;

      std::string name           = j->first,
        units               = diag->get_metadata().get_string("units"),
        glaciological_units = diag->get_metadata().get_string("glaciological_units");

      if (glaciological_units.empty() == false) {
        units = glaciological_units;
      }

      if (diag->get_metadata().get_n_spatial_dimensions() == d) {

        PetscPrintf(grid.com, "   Name: %s [%s]\n", name.c_str(), units.c_str());

        for (int k = 0; k < diag->get_nvars(); ++k) {
          NCSpatialVariable var = diag->get_metadata(k);

          std::string long_name = var.get_string("long_name");

          PetscPrintf(grid.com, "      -  %s\n", long_name.c_str());
        }

        PetscPrintf(grid.com, "\n");

      }

      ++j;
    }
  }

  // scalar time-series
  PetscPrintf(grid.com, "======== Available time-series ========\n");

  std::map<std::string, TSDiagnostic*>::iterator j = ts_diagnostics.begin();
  while (j != ts_diagnostics.end()) {
    TSDiagnostic *diag = j->second;

    std::string name = j->first,
      long_name = diag->get_string("long_name"),
      units = diag->get_string("units"),
      glaciological_units = diag->get_string("glaciological_units");

    if (glaciological_units.empty() == false) {
      units = glaciological_units;
    }

    PetscPrintf(grid.com,
                "   Name: %s [%s]\n"
                "      -  %s\n\n",
                name.c_str(), units.c_str(), long_name.c_str());

    ++j;
  }
}


IceModel_hardav::IceModel_hardav(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "hardav", grid));

  // choice to use SSA power; see #285
  const double power = 1.0 / grid.config.get("ssa_Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  vars[0].set_double("valid_min", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

//! \brief Computes vertically-averaged ice hardness.
void IceModel_hardav::compute(IceModelVec* &output) {
  const double fillval = grid.config.get("fill_value");
  double *Eij; // columns of enthalpy values

  const IceFlowLaw *flow_law = model->stress_balance->get_stressbalance()->get_flow_law();
  if (flow_law == NULL) {
    flow_law = model->stress_balance->get_ssb_modifier()->get_flow_law();
    if (flow_law == NULL) {
      throw RuntimeError("Can't compute vertically-averaged hardness: no flow law is used.");
    }
  }

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "hardav", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(model->Enth3);
  list.add(model->ice_thickness);
  list.add(*result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    model->Enth3.getInternalColumn(i,j,&Eij);
    const double H = model->ice_thickness(i,j);
    if (mask.icy(i, j)) {
      (*result)(i,j) = flow_law->averaged_hardness(H, grid.kBelowHeight(H),
                                                   &(grid.z()[0]), Eij);
    } else { // put negative value below valid range
      (*result)(i,j) = fillval;
    }
  }

  output = result;
}


IceModel_rank::IceModel_rank(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "rank", grid));

  set_attrs("processor rank", "", "1", "", 0);
  vars[0].set_time_independent(true);
}

void IceModel_rank::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "rank", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  IceModelVec::AccessList list;
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    (*result)(p.i(),p.j()) = grid.rank();
  }

  output = result;
}


IceModel_cts::IceModel_cts(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "cts", grid, g.z()));

  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

void IceModel_cts::compute(IceModelVec* &output) {

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.z());

  IceModelVec3 *result = new IceModelVec3;
  result->create(grid, "cts", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  model->setCTSFromEnthalpy(*result);

  output = result;
}

IceModel_proc_ice_area::IceModel_proc_ice_area(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "proc_ice_area", grid));

  set_attrs("number of cells containing ice in a processor's domain", "",
            "", "", 0);
  vars[0].set_time_independent(true);
}

void IceModel_proc_ice_area::compute(IceModelVec* &output) {

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");
  IceModelVec2Int *ice_mask = variables.get_2d_mask("mask");

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "proc_ice_area", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  int ice_filled_cells = 0;

  MaskQuery mask(*ice_mask);

  IceModelVec::AccessList list;
  list.add(*ice_mask);
  list.add(*thickness);
  for (Points p(grid); p; p.next()) {
    if (mask.icy(p.i(), p.j())) {
      ice_filled_cells += 1;
    }
  }

  for (Points p(grid); p; p.next()) {
    (*result)(p.i(), p.j()) = ice_filled_cells;
  }

  output = result;
}


IceModel_temp::IceModel_temp(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "temp", grid, g.z()));

  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  vars[0].set_double("valid_min", 0);
}

void IceModel_temp::compute(IceModelVec* &output) {

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.z());

  IceModelVec3 *result = new IceModelVec3;
  result->create(grid, "temp", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");
  IceModelVec3 *enthalpy = variables.get_3d_scalar("enthalpy");

  double *Tij, *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result->getInternalColumn(i,j,&Tij);
    enthalpy->getInternalColumn(i,j,&Enthij);
    for (unsigned int k=0; k <grid.Mz(); ++k) {
      const double depth = (*thickness)(i,j) - grid.z(k);
      Tij[k] = model->EC->getAbsTemp(Enthij[k],
                                     model->EC->getPressureFromDepth(depth));
    }
  }

  output = result;
}


IceModel_temp_pa::IceModel_temp_pa(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "temp_pa", grid, g.z()));

  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
  vars[0].set_double("valid_max", 0);
}

void IceModel_temp_pa::compute(IceModelVec* &output) {
  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");
  double melting_point_temp = grid.config.get("water_melting_point_temperature");

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.z());

  IceModelVec3 *result = new IceModelVec3;
  result->create(grid, "temp_pa", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");
  IceModelVec3  *enthalpy  = variables.get_3d_scalar("enthalpy");

  double *Tij, *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  for (Points pt(grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    result->getInternalColumn(i,j,&Tij);
    enthalpy->getInternalColumn(i,j,&Enthij);
    for (unsigned int k=0; k < grid.Mz(); ++k) {
      const double depth = (*thickness)(i,j) - grid.z(k),
        p = model->EC->getPressureFromDepth(depth);
      Tij[k] = model->EC->getPATemp(Enthij[k], p);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
        // is 273.15
        if (model->EC->isTemperate(Enthij[k],p) && ((*thickness)(i,j) > 0)) {
          Tij[k] = melting_point_temp;
        }
      }

    }
  }

  result->shift(-melting_point_temp);

  output = result;
}

IceModel_temppabase::IceModel_temppabase(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "temppabase", grid));

  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "Celsius", "Celsius", 0);
}

void IceModel_temppabase::compute(IceModelVec* &output) {

  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");
  double melting_point_temp = grid.config.get("water_melting_point_temperature");

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "temp_pa_base", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");
  IceModelVec3 *enthalpy = variables.get_3d_scalar("enthalpy");

  double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  for (Points pt(grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    enthalpy->getInternalColumn(i,j,&Enthij);

    const double depth = (*thickness)(i,j),
      p = model->EC->getPressureFromDepth(depth);
    (*result)(i,j) = model->EC->getPATemp(Enthij[0], p);

    if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
      // is 273.15
      if (model->EC->isTemperate(Enthij[0],p) && ((*thickness)(i,j) > 0)) {
        (*result)(i,j) = melting_point_temp;
      }
    }
  }

  result->shift(-melting_point_temp);

  output = result;
}

IceModel_enthalpysurf::IceModel_enthalpysurf(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "enthalpysurf", grid));

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

void IceModel_enthalpysurf::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "enthalpysurf", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  double fill_value = grid.config.get("fill_value");

  // compute levels corresponding to 1 m below the ice surface:

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) = std::max(model->ice_thickness(i,j) - 1.0, 0.0);
  }

  model->Enth3.getSurfaceValues(*result, *result);  // z=0 slice

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (model->ice_thickness(i,j) <= 1.0) {
      (*result)(i,j) = fill_value;
    }
  }

  output = result;
}

IceModel_enthalpybase::IceModel_enthalpybase(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "enthalpybase", grid));

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

void IceModel_enthalpybase::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "enthalpybase", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  model->Enth3.getHorSlice(*result, 0.0);  // z=0 slice

  result->mask_by(model->ice_thickness, grid.config.get("fill_value"));

  output = result;
}


IceModel_tempbase::IceModel_tempbase(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "tempbase", grid));

  set_attrs("ice temperature at the base of ice", "",
            "K", "K", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

void IceModel_tempbase::compute(IceModelVec* &output) {

  IceModelVec2S *result = NULL,
    *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModel_enthalpybase enth(model, grid, variables);

  enth.compute(output);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  // result contains basal enthalpy; note that it is allocated by
  // enth.compute().

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(*thickness);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double depth = (*thickness)(i,j),
      pressure = model->EC->getPressureFromDepth(depth);
    if (mask.icy(i, j)) {
      (*result)(i,j) = model->EC->getAbsTemp((*result)(i,j), pressure);
    } else {
      (*result)(i,j) = grid.config.get("fill_value");
    }
  }


  result->metadata() = vars[0];
  output = result;
}

IceModel_tempsurf::IceModel_tempsurf(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "tempsurf", grid));

  set_attrs("ice temperature at 1m below the ice surface", "",
            "K", "K", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

void IceModel_tempsurf::compute(IceModelVec* &output) {

  IceModelVec2S *result, *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModel_enthalpysurf enth(model, grid, variables);

  enth.compute(output);
  result = dynamic_cast<IceModelVec2S*>(output);
  if (result == NULL) {
    throw RuntimeError( "dynamic_cast failure");
  }

  // result contains surface enthalpy; note that it is allocated by
  // enth.compute().

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  double depth = 1.0,
    pressure = model->EC->getPressureFromDepth(depth);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*thickness)(i,j) > 1) {
      (*result)(i,j) = model->EC->getAbsTemp((*result)(i,j), pressure);
    } else {
      (*result)(i,j) = grid.config.get("fill_value");
    }
  }


  result->metadata() = vars[0];
  output = result;
}


IceModel_liqfrac::IceModel_liqfrac(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "liqfrac", grid, g.z()));

  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
  vars[0].set_double("valid_min", 0);
  vars[0].set_double("valid_max", 1);
}

void IceModel_liqfrac::compute(IceModelVec* &output) {

  // update vertical levels (in case the grid was extended
  vars[0].set_levels(grid.z());

  IceModelVec3 *result = new IceModelVec3;
  result->create(grid, "liqfrac", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  bool cold_mode = grid.config.get_flag("do_cold_ice_methods");

  if (cold_mode) {
    result->set(0.0);
  } else {
    model->compute_liquid_water_fraction(model->Enth3, *result);
  }

  output = result;
}

IceModel_tempicethk::IceModel_tempicethk(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "tempicethk", grid));

  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

void IceModel_tempicethk::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "tempicethk", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  double *Enth;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(model->Enth3);
  list.add(model->ice_thickness);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      model->Enth3.getInternalColumn(i,j,&Enth);
      double temperate_ice_thickness = 0.0;
      double ice_thickness = model->ice_thickness(i,j);
      const unsigned int ks = grid.kBelowHeight(ice_thickness);

      for (unsigned int k=0; k<ks; ++k) { // FIXME issue #15
        double pressure = model->EC->getPressureFromDepth(ice_thickness - grid.z(k));

        if (model->EC->isTemperate(Enth[k], pressure)) {
          temperate_ice_thickness += grid.z(k+1) - grid.z(k);
        }
      }

      double pressure = model->EC->getPressureFromDepth(ice_thickness - grid.z(ks));
      if (model->EC->isTemperate(Enth[ks], pressure)) {
        temperate_ice_thickness += ice_thickness - grid.z(ks);
      }

      (*result)(i,j) = temperate_ice_thickness;
    } else {
      // ice-free
      (*result)(i,j) = grid.config.get("fill_value");
    }
  }

  output = result;
}

IceModel_tempicethk_basal::IceModel_tempicethk_basal(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "tempicethk_basal", grid));

  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
void IceModel_tempicethk_basal::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "tempicethk_basal", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  double *Enth;
  EnthalpyConverter *EC = model->EC;

  MaskQuery mask(model->vMask);

  const double fill_value = grid.config.get("fill_value");

  IceModelVec::AccessList list;
  list.add(model->vMask);
  list.add(*result);
  list.add(model->ice_thickness);
  list.add(model->Enth3);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double thk = model->ice_thickness(i,j);

    // if we have no ice, go on to the next grid point (this cell will be
    // marked as "missing" later)
    if (mask.ice_free(i, j)) {
      (*result)(i,j) = fill_value;
      continue;
    }

    model->Enth3.getInternalColumn(i,j,&Enth);
    double pressure;
    unsigned int ks = grid.kBelowHeight(thk),
      k = 0;

    while (k <= ks) {         // FIXME issue #15
      pressure = EC->getPressureFromDepth(thk - grid.z(k));

      if (EC->isTemperate(Enth[k],pressure)) {
        k++;
      } else {
        break;
      }
    }
    // after this loop 'pressure' is equal to the pressure at the first level
    // that is cold

    // no temperate ice at all; go to the next grid point
    if (k == 0) {
      (*result)(i,j) = 0.0;
      continue;
    }

    // the whole column is temperate (except, possibly, some ice between
    // z(ks) and the total thickness; we ignore it)
    if (k == ks + 1) {
      (*result)(i,j) = grid.z(ks);
      continue;
    }

    double
      pressure_0 = EC->getPressureFromDepth(thk - grid.z(k-1)),
      dz         = grid.z(k) - grid.z(k-1),
      slope1     = (Enth[k] - Enth[k-1]) / dz,
      slope2     = (EC->getEnthalpyCTS(pressure) - EC->getEnthalpyCTS(pressure_0)) / dz;

    if (slope1 != slope2) {
      (*result)(i,j) = grid.z(k-1) +
        (EC->getEnthalpyCTS(pressure_0) - Enth[k-1]) / (slope1 - slope2);

      // check if the resulting thickness is valid:
      (*result)(i,j) = std::max((*result)(i,j), grid.z(k-1));
      (*result)(i,j) = std::min((*result)(i,j), grid.z(k));
    } else {
      throw RuntimeError::formatted("Linear interpolation of the thickness of the basal temperate layer failed:\n"
                                    "(i=%d, j=%d, k=%d, ks=%d)\n",
                                    i, j, k, ks);
    }
  }

  output = result;
}

IceModel_flux_divergence::IceModel_flux_divergence(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "flux_divergence", grid));

  set_attrs("flux divergence", "", "m s-1", "m year-1", 0);
}

void IceModel_flux_divergence::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "flux_divergence", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->flux_divergence);

  output = result;
}

IceModel_climatic_mass_balance_cumulative::IceModel_climatic_mass_balance_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "climatic_mass_balance_cumulative", grid));

  set_attrs("cumulative ice-equivalent climatic mass balance", "",
            "kg m-2", "kg m-2", 0);
}

void IceModel_climatic_mass_balance_cumulative::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "climatic_mass_balance_cumulative", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->climatic_mass_balance_cumulative);

  output = result;
}

IceModel_ivol::IceModel_ivol(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivol", time_dimension_name);

  ts->get_metadata().set_units("m3");
  ts->get_dimension_metadata().set_units(time_units);

  ts->get_metadata().set_string("long_name", "total ice volume");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_ivol::update(double a, double b) {
  double value;

  model->compute_ice_volume(value);

  ts->append(value, a, b);
}

IceModel_slvol::IceModel_slvol(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "slvol", time_dimension_name);

  ts->get_metadata().set_units("m");
  ts->get_dimension_metadata().set_units(time_units);

  ts->get_metadata().set_string("long_name", "total sea-level relevant ice IN SEA-LEVEL EQUIVALENT");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_slvol::update(double a, double b) {
  double value;

  model->compute_sealevel_volume(value);

  ts->append(value, a, b);
}

IceModel_divoldt::IceModel_divoldt(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "divoldt", time_dimension_name);

  ts->get_metadata().set_units("m3 s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->rate_of_change = true;

  ts->get_metadata().set_string("long_name", "total ice volume rate of change");
}

void IceModel_divoldt::update(double a, double b) {
  double value;

  model->compute_ice_volume(value);

  // note that "value" below *should* be the ice volume
  ts->append(value, a, b);
}


IceModel_iarea::IceModel_iarea(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iarea", time_dimension_name);

  ts->get_metadata().set_units("m2");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total ice area");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_iarea::update(double a, double b) {
  double value;

  model->compute_ice_area(value);

  ts->append(value, a, b);
}

IceModel_imass::IceModel_imass(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "imass", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total ice mass");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_imass::update(double a, double b) {
  double value;

  model->compute_ice_volume(value);

  ts->append(value * grid.config.get("ice_density"), a, b);
}


IceModel_dimassdt::IceModel_dimassdt(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "dimassdt", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total ice mass rate of change");

  ts->rate_of_change = true;
}

void IceModel_dimassdt::update(double a, double b) {
  double value;

  model->compute_ice_volume(value);

  ts->append(value * grid.config.get("ice_density"), a, b);
}


IceModel_ivoltemp::IceModel_ivoltemp(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivoltemp", time_dimension_name);

  ts->get_metadata().set_units("m3");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total volume of temperate ice");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_ivoltemp::update(double a, double b) {
  double value;

  model->compute_ice_volume_temperate(value);

  ts->append(value, a, b);
}


IceModel_ivolcold::IceModel_ivolcold(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolcold", time_dimension_name);

  ts->get_metadata().set_units("m3");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total volume of cold ice");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_ivolcold::update(double a, double b) {
  double value;

  model->compute_ice_volume_cold(value);

  ts->append(value, a, b);
}

IceModel_iareatemp::IceModel_iareatemp(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareatemp", time_dimension_name);

  ts->get_metadata().set_units("m2");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "ice-covered area where basal ice is temperate");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_iareatemp::update(double a, double b) {
  double value;

  model->compute_ice_area_temperate(value);

  ts->append(value, a, b);
}

IceModel_iareacold::IceModel_iareacold(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareacold", time_dimension_name);

  ts->get_metadata().set_units("m2");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "ice-covered area where basal ice is cold");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_iareacold::update(double a, double b) {
  double value;

  model->compute_ice_area_cold(value);

  ts->append(value, a, b);
}

IceModel_ienthalpy::IceModel_ienthalpy(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ienthalpy", time_dimension_name);

  ts->get_metadata().set_units("J");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total ice enthalpy");
  ts->get_metadata().set_double("valid_min", 0.0);
}

void IceModel_ienthalpy::update(double a, double b) {
  double value;

  model->compute_ice_enthalpy(value);

  ts->append(value, a, b);
}

IceModel_iareag::IceModel_iareag(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareag", time_dimension_name);

  ts->get_metadata().set_units("m2");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total grounded ice area");
}

void IceModel_iareag::update(double a, double b) {
  double value;

  model->compute_ice_area_grounded(value);

  ts->append(value, a, b);
}

IceModel_iareaf::IceModel_iareaf(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "iareaf", time_dimension_name);

  ts->get_metadata().set_units("m2");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total floating ice area");
}

void IceModel_iareaf::update(double a, double b) {
  double value;

  model->compute_ice_area_floating(value);

  ts->append(value, a, b);
}

IceModel_dt::IceModel_dt(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "dt", time_dimension_name);

  ts->get_metadata().set_units("second");
  ts->get_metadata().set_glaciological_units("year");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "mass continuity time step");
}

void IceModel_dt::update(double a, double b) {

  ts->append(model->dt, a, b);
}

IceModel_max_diffusivity::IceModel_max_diffusivity(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "max_diffusivity", time_dimension_name);

  ts->get_metadata().set_units("m2 s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "maximum diffusivity");
}

void IceModel_max_diffusivity::update(double a, double b) {
  double value;

  model->stress_balance->get_max_diffusivity(value);

  ts->append(value, a, b);
}

IceModel_surface_flux::IceModel_surface_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "surface_ice_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total over ice domain of top surface ice mass flux");
  ts->rate_of_change = true;
}

void IceModel_surface_flux::update(double a, double b) {
  double value;

  value = model->surface_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_surface_flux_cumulative::IceModel_surface_flux_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "surface_ice_flux_cumulative", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "cumulative total over ice domain of top surface ice mass flux");
}

void IceModel_surface_flux_cumulative::update(double a, double b) {
  double value;

  value = model->surface_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_grounded_basal_flux::IceModel_grounded_basal_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "grounded_basal_ice_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total over grounded ice domain of basal mass flux");
  ts->rate_of_change = true;
}

void IceModel_grounded_basal_flux::update(double a, double b) {
  double value;

  value = model->grounded_basal_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_grounded_basal_flux_cumulative::IceModel_grounded_basal_flux_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "grounded_basal_ice_flux_cumulative", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "cumulative total grounded basal mass flux");
}

void IceModel_grounded_basal_flux_cumulative::update(double a, double b) {
  double value;

  value = model->grounded_basal_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_sub_shelf_flux::IceModel_sub_shelf_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sub_shelf_ice_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total sub-shelf ice flux");
  ts->rate_of_change = true;
}

void IceModel_sub_shelf_flux::update(double a, double b) {
  double value;

  value = model->sub_shelf_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_sub_shelf_flux_cumulative::IceModel_sub_shelf_flux_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sub_shelf_ice_flux_cumulative", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "cumulative total sub-shelf ice flux");
}

void IceModel_sub_shelf_flux_cumulative::update(double a, double b) {
  double value;

  value = model->sub_shelf_ice_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_nonneg_flux::IceModel_nonneg_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "nonneg_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  ts->rate_of_change = true;
}

void IceModel_nonneg_flux::update(double a, double b) {
  double value;

  value = model->nonneg_rule_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_nonneg_flux_cumulative::IceModel_nonneg_flux_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "nonneg_flux_cumulative", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "cumulative 'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
}

void IceModel_nonneg_flux_cumulative::update(double a, double b) {
  double value;

  value = model->nonneg_rule_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_discharge_flux::IceModel_discharge_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "discharge_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "discharge (calving & icebergs) flux");
  ts->rate_of_change = true;
}

void IceModel_discharge_flux::update(double a, double b) {
  double value = model->discharge_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_discharge_flux_cumulative::IceModel_discharge_flux_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "discharge_flux_cumulative", time_dimension_name);

  ts->get_metadata().set_units("kg");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "cumulative discharge (calving etc.) flux");
}

void IceModel_discharge_flux_cumulative::update(double a, double b) {
  double value = model->discharge_flux_cumulative;

  ts->append(value, a, b);
}

IceModel_dHdt::IceModel_dHdt(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "dHdt", grid));

  set_attrs("ice thickness rate of change", "tendency_of_land_ice_thickness",
            "m s-1", "m year-1", 0);

  vars[0].set_double("valid_min",  grid.convert(-1e6, "m/year", "m/s"));
  vars[0].set_double("valid_max",  grid.convert( 1e6, "m/year", "m/s"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
  vars[0].set_string("cell_methods", "time: mean");

  last_ice_thickness.create(grid, "last_ice_thickness", WITHOUT_GHOSTS);
  last_ice_thickness.set_attrs("internal",
                               "ice thickness at the time of the last report of dHdt",
                               "m", "land_ice_thickness");

  last_report_time = GSL_NAN;
}

void IceModel_dHdt::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "dHdt", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  if (gsl_isnan(last_report_time)) {
    result->set(grid.convert(2e6, "m/year", "m/s"));
  } else {
    IceModelVec::AccessList list;
    list.add(*result);
    list.add(last_ice_thickness);
    list.add(model->ice_thickness);

    double dt = grid.time->current() - last_report_time;
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (model->ice_thickness(i, j) - last_ice_thickness(i, j)) / dt;
    }
  }

  // Save the ice thickness and the corresponding time:
  this->update_cumulative();

  output = result;
}

void IceModel_dHdt::update_cumulative() {
  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(last_ice_thickness);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    last_ice_thickness(i, j) = model->ice_thickness(i, j);
  }

  last_report_time = grid.time->current();
}

IceModel_ivolg::IceModel_ivolg(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolg", time_dimension_name);

  ts->get_metadata().set_units("m3");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total grounded ice volume");
}

void IceModel_ivolg::update(double a, double b) {
  double volume=0.0, value;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(model->vMask);
  list.add(model->cell_area);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.grounded_ice(i,j)) {
      volume += model->cell_area(i,j) * model->ice_thickness(i,j);
    }
  }

  GlobalSum(grid.com, &volume,  &value);

  ts->append(value, a, b);
}

IceModel_ivolf::IceModel_ivolf(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "ivolf", time_dimension_name);

  ts->get_metadata().set_units("m3");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "total floating ice volume");
}

void IceModel_ivolf::update(double a, double b) {
  double volume=0.0, value;

  MaskQuery mask(model->vMask);

  IceModelVec::AccessList list;
  list.add(model->ice_thickness);
  list.add(model->vMask);
  list.add(model->cell_area);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.floating_ice(i,j)) {
      volume += model->cell_area(i,j) * model->ice_thickness(i,j);
    }
  }

  GlobalSum(grid.com, &volume,  &value);

  ts->append(value, a, b);
}

//! \brief Reports the maximum horizontal absolute velocity component over the grid.
/*!
 * This is the value used by the adaptive time-stepping code in the CFL condition
 * for horizontal advection (i.e. in energy and mass conservation time steps).
 *
 * This is not the maximum horizontal speed, but rather the maximum of components.
 *
 * Note that this picks up the value computed during the time-step taken at a
 * reporting time. (It is not the "average over the reporting interval computed using
 * differencing in time", as other rate-of-change diagnostics.)
 */
IceModel_max_hor_vel::IceModel_max_hor_vel(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {

  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "max_hor_vel", time_dimension_name);

  ts->get_metadata().set_units("m/second");
  ts->get_metadata().set_glaciological_units("m/year");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name",
                                "maximum abs component of horizontal ice velocity"
                                " over grid in last time step during time-series reporting interval");
}

void IceModel_max_hor_vel::update(double a, double b) {

  double gmaxu = model->gmaxu, gmaxv = model->gmaxv;

  ts->append(gmaxu > gmaxv ? gmaxu : gmaxv, a, b);
}

IceModel_H_to_Href_flux::IceModel_H_to_Href_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "H_to_Href_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "mass flux from thk to Href");
  ts->rate_of_change = true;
}

void IceModel_H_to_Href_flux::update(double a, double b) {

  ts->append(model->H_to_Href_flux_cumulative, a, b);
}


IceModel_Href_to_H_flux::IceModel_Href_to_H_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "Href_to_H_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "mass flux from Href to thk");
  ts->rate_of_change = true;
}

void IceModel_Href_to_H_flux::update(double a, double b) {

  ts->append(model->Href_to_H_flux_cumulative, a, b);
}



IceModel_sum_divQ_flux::IceModel_sum_divQ_flux(IceModel *m, IceGrid &g, Vars &my_vars)
  : TSDiag<IceModel>(m, g, my_vars) {
  // set metadata:
  ts = new DiagnosticTimeseries(&grid, "sum_divQ_flux", time_dimension_name);

  ts->get_metadata().set_units("kg s-1");
  ts->get_dimension_metadata().set_units(time_units);
  ts->get_metadata().set_string("long_name", "sum(divQ)");
  ts->rate_of_change = true;
}

void IceModel_sum_divQ_flux::update(double a, double b) {

  ts->append(model->sum_divQ_SIA_cumulative + model->sum_divQ_SSA_cumulative,
             a, b);
}


IceModel_nonneg_flux_2D_cumulative::IceModel_nonneg_flux_2D_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "nonneg_flux_cumulative", grid));

  set_attrs("cumulative non-negative rule (thk >= 0) flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

void IceModel_nonneg_flux_2D_cumulative::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->nonneg_flux_2D_cumulative);

  output = result;
}

IceModel_grounded_basal_flux_2D_cumulative::IceModel_grounded_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "grounded_basal_flux_cumulative", grid));

  set_attrs("cumulative grounded basal flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

void IceModel_grounded_basal_flux_2D_cumulative::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->grounded_basal_flux_2D_cumulative);

  output = result;
}

IceModel_floating_basal_flux_2D_cumulative::IceModel_floating_basal_flux_2D_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "floating_basal_flux_cumulative", grid));

  set_attrs("cumulative floating basal flux (positive means ice gain)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

void IceModel_floating_basal_flux_2D_cumulative::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->floating_basal_flux_2D_cumulative);

  output = result;
}


IceModel_discharge_flux_2D_cumulative::IceModel_discharge_flux_2D_cumulative(IceModel *m, IceGrid &g, Vars &my_vars)
  : Diag<IceModel>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), "discharge_flux_cumulative", grid));

  set_attrs("cumulative ice discharge (calving) flux (negative means ice loss)",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
}

void IceModel_discharge_flux_2D_cumulative::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->discharge_flux_2D_cumulative);

  output = result;
}

#if (PISM_USE_PROJ4==1)
IceModel_lat_lon_bounds::IceModel_lat_lon_bounds(IceModel *m, IceGrid &g, Vars &my_vars,
                                                 std::string var_name, std::string proj_string)
  : Diag<IceModel>(m, g, my_vars) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  std::vector<double> levels(4);
  for (int k = 0; k < 4; ++k) {
    levels[k] = k;
  }

  vars.push_back(NCSpatialVariable(grid.config.get_unit_system(), m_var_name + "_bnds", grid, levels));
  vars[0].get_z().set_name("nv4");
  vars[0].get_z().clear_all_strings();
  vars[0].get_z().clear_all_doubles();
  vars[0].set_time_independent(true);

  if (m_var_name == "lon") {
    set_attrs("longitude bounds", "", "degree_east", "degree_east", 0);
    vars[0].set_double("valid_min", -180);
    vars[0].set_double("valid_max", 180);
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    vars[0].set_double("valid_min", -90);
    vars[0].set_double("valid_max", 90);
  }

  lonlat = pj_init_plus("+proj=latlong +datum=WGS84 +ellps=WGS84");
  if (lonlat == NULL) {
    throw RuntimeError("projection initialization failed\n"
                       "('+proj=latlong +datum=WGS84 +ellps=WGS84').\n");
  }

  pism = pj_init_plus(proj_string.c_str());
  if (pism == NULL) {
    throw RuntimeError::formatted("proj.4 string '%s' is invalid.", proj_string.c_str());
  }
}

IceModel_lat_lon_bounds::~IceModel_lat_lon_bounds() {
  pj_free(pism);
  pj_free(lonlat);
}

void IceModel_lat_lon_bounds::compute(IceModelVec* &output) {

  IceModelVec3Custom *result = new IceModelVec3Custom;
  std::map<std::string,std::string> attrs;
  std::vector<double> indices(4);

  result->create(grid, m_var_name + "_bnds", "nv4",
                 indices, attrs);
  result->metadata() = vars[0];

  double dx2 = 0.5 * grid.dx(), dy2 = 0.5 * grid.dy();
  double x_offsets[] = {-dx2, dx2, dx2, -dx2};
  double y_offsets[] = {-dy2, -dy2, dy2, dy2};

  bool latitude = true;
  if (m_var_name == "lon") {
    latitude = false;
  }

  IceModelVec::AccessList list;
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double x0 = grid.x(i), y0 = grid.y(j);
    double *values;

    result->getInternalColumn(i,j,&values);

    for (int k = 0; k < 4; ++k) {
      double
        x = x0 + x_offsets[k],
        y = y0 + y_offsets[k];

      // compute lon,lat coordinates:
      pj_transform(pism, lonlat, 1, 1, &x, &y, NULL);

      // NB! proj.4 converts x,y pairs into lon,lat pairs in *radians*.

      if (latitude) {
        values[k] = y * RAD_TO_DEG;
      } else {
        values[k] = x * RAD_TO_DEG;
      }
    }
  }

  output = result;
}
#elif (PISM_USE_PROJ4==0)
  // do nothing
#else  // PISM_USE_PROJ4 is not set
#error "PISM build system error: PISM_USE_PROJ4 is not set."
#endif

} // end of namespace pism
