// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include <gsl/gsl_math.h>

#include "iceModel_diagnostics.hh"

#include "base/rheology/FlowLaw.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/SSB_Modifier.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec3Custom.hh"
#include "enthalpyConverter.hh"
#include "base/util/PISMVars.hh"
#include "base/util/pism_utilities.hh"
#include "coupler/PISMOcean.hh"
#include "earth/PISMbedDef.hh"

#include "base/grounded_cell_fraction.hh"
#include "base/part_grid_threshold_thickness.hh"
#include "base/util/projection.hh"
#include "base/energy/utilities.hh"
#include "base/energy/EnergyModel.hh"

#if (PISM_USE_PROJ4==1)
#include "base/util/Proj.hh"
#endif

namespace pism {

// Horrendous names used by InitMIP (and ISMIP6, and CMIP5). Ugh.
static const char* land_ice_area_fraction_name           = "sftgif";
static const char* grounded_ice_sheet_area_fraction_name = "sftgrf";
static const char* floating_ice_sheet_area_fraction_name = "sftflf";

void IceModel::init_diagnostics() {

  // Add IceModel diagnostics:
  m_diagnostics["cts"]              = Diagnostic::Ptr(new IceModel_cts(this));
  m_diagnostics["enthalpybase"]     = Diagnostic::Ptr(new IceModel_enthalpybase(this));
  m_diagnostics["enthalpysurf"]     = Diagnostic::Ptr(new IceModel_enthalpysurf(this));
  m_diagnostics["hardav"]           = Diagnostic::Ptr(new IceModel_hardav(this));
  m_diagnostics["hardness"]         = Diagnostic::Ptr(new IceModel_hardness(this));
  m_diagnostics["liqfrac"]          = Diagnostic::Ptr(new IceModel_liqfrac(this));
  m_diagnostics["proc_ice_area"]    = Diagnostic::Ptr(new IceModel_proc_ice_area(this));
  m_diagnostics["rank"]             = Diagnostic::Ptr(new IceModel_rank(this));
  m_diagnostics["temp"]             = Diagnostic::Ptr(new IceModel_temp(this));
  m_diagnostics["temp_pa"]          = Diagnostic::Ptr(new IceModel_temp_pa(this));
  m_diagnostics["tempbase"]         = Diagnostic::Ptr(new IceModel_tempbase(this));
  m_diagnostics["tempicethk"]       = Diagnostic::Ptr(new IceModel_tempicethk(this));
  m_diagnostics["tempicethk_basal"] = Diagnostic::Ptr(new IceModel_tempicethk_basal(this));
  m_diagnostics["temppabase"]       = Diagnostic::Ptr(new IceModel_temppabase(this));
  m_diagnostics["tempsurf"]         = Diagnostic::Ptr(new IceModel_tempsurf(this));
  m_diagnostics["dHdt"]             = Diagnostic::Ptr(new IceModel_dHdt(this));
  m_diagnostics["effective_viscosity"] = Diagnostic::Ptr(new IceModel_viscosity(this));

  m_diagnostics[land_ice_area_fraction_name]           = Diagnostic::Ptr(new IceModel_land_ice_area_fraction(this));
  m_diagnostics[grounded_ice_sheet_area_fraction_name] = Diagnostic::Ptr(new IceModel_grounded_ice_sheet_area_fraction(this));
  m_diagnostics[floating_ice_sheet_area_fraction_name] = Diagnostic::Ptr(new IceModel_floating_ice_sheet_area_fraction(this));

  m_diagnostics["flux_divergence"]                  = Diagnostic::Ptr(new IceModel_flux_divergence(this));
  m_diagnostics["climatic_mass_balance_cumulative"] = Diagnostic::Ptr(new IceModel_climatic_mass_balance_cumulative(this));
  m_diagnostics["nonneg_flux_cumulative"]           = Diagnostic::Ptr(new IceModel_nonneg_flux_2D_cumulative(this));
  m_diagnostics["grounded_basal_flux_cumulative"]   = Diagnostic::Ptr(new IceModel_grounded_basal_flux_2D_cumulative(this));
  m_diagnostics["floating_basal_flux_cumulative"]   = Diagnostic::Ptr(new IceModel_floating_basal_flux_2D_cumulative(this));
  m_diagnostics["discharge_flux_cumulative"]        = Diagnostic::Ptr(new IceModel_discharge_flux_2D_cumulative(this));
  m_diagnostics["discharge_flux"]                   = Diagnostic::Ptr(new IceModel_discharge_flux_2D(this));

  m_diagnostics["surface_mass_balance_average"] = Diagnostic::Ptr(new IceModel_surface_mass_balance_average(this));
  m_diagnostics["basal_mass_balance_average"]   = Diagnostic::Ptr(new IceModel_basal_mass_balance_average(this));
  m_diagnostics["height_above_flotation"]       = Diagnostic::Ptr(new IceModel_height_above_flotation(this));
  m_diagnostics["ice_mass"]                     = Diagnostic::Ptr(new IceModel_ice_mass(this));
  m_diagnostics["topg_sl_adjusted"]             = Diagnostic::Ptr(new IceModel_topg_sl_adjusted(this));

#if (PISM_USE_PROJ4==1)
  std::string proj4 = m_grid->get_mapping_info().proj4;
  if (not proj4.empty()) {
    m_diagnostics["lat_bnds"] = Diagnostic::Ptr(new IceModel_lat_lon_bounds(this, "lat", proj4));
    m_diagnostics["lon_bnds"] = Diagnostic::Ptr(new IceModel_lat_lon_bounds(this, "lon", proj4));
  }
#endif

  m_ts_diagnostics["volume_glacierized"]                   = TSDiagnostic::Ptr(new IceModel_volume_glacierized(this));
  m_ts_diagnostics["volume_nonglacierized"]                = TSDiagnostic::Ptr(new IceModel_volume_nonglacierized(this));
  m_ts_diagnostics["slvol"]                                = TSDiagnostic::Ptr(new IceModel_slvol(this));
  m_ts_diagnostics["volume_rate_of_change_glacierized"]    = TSDiagnostic::Ptr(new IceModel_volume_rate_of_change_glacierized(this));
  m_ts_diagnostics["volume_rate_of_change_nonglacierized"] = TSDiagnostic::Ptr(new IceModel_volume_rate_of_change_nonglacierized(this));
  m_ts_diagnostics["area_glacierized"]                     = TSDiagnostic::Ptr(new IceModel_area_glacierized(this));
  m_ts_diagnostics["mass_glacierized"]                     = TSDiagnostic::Ptr(new IceModel_mass_glacierized(this));
  m_ts_diagnostics["mass_nonglacierized"]                  = TSDiagnostic::Ptr(new IceModel_mass_nonglacierized(this));
  m_ts_diagnostics["mass_rate_of_change_glacierized"]      = TSDiagnostic::Ptr(new IceModel_mass_rate_of_change_glacierized(this));
  m_ts_diagnostics["mass_rate_of_change_nonglacierized"]   = TSDiagnostic::Ptr(new IceModel_mass_rate_of_change_nonglacierized(this));
  m_ts_diagnostics["volume_glacierized_temperate"]         = TSDiagnostic::Ptr(new IceModel_volume_glacierized_temperate(this));
  m_ts_diagnostics["volume_nonglacierized_temperate"]      = TSDiagnostic::Ptr(new IceModel_volume_nonglacierized_temperate(this));
  m_ts_diagnostics["volume_glacierized_cold"]              = TSDiagnostic::Ptr(new IceModel_volume_glacierized_cold(this));
  m_ts_diagnostics["volume_nonglacierized_cold"]           = TSDiagnostic::Ptr(new IceModel_volume_nonglacierized_cold(this));
  m_ts_diagnostics["volume_glacierized_grounded"]          = TSDiagnostic::Ptr(new IceModel_volume_glacierized_grounded(this));
  m_ts_diagnostics["volume_glacierized_shelf"]             = TSDiagnostic::Ptr(new IceModel_volume_glacierized_shelf(this));
  m_ts_diagnostics["area_glacierized_temperate_base"]      = TSDiagnostic::Ptr(new IceModel_area_glacierized_temperate_base(this));
  m_ts_diagnostics["area_glacierized_cold_base"]           = TSDiagnostic::Ptr(new IceModel_area_glacierized_cold_base(this));
  m_ts_diagnostics["area_glacierized_grounded"]            = TSDiagnostic::Ptr(new IceModel_area_glacierized_grounded(this));
  m_ts_diagnostics["area_glacierized_shelf"]               = TSDiagnostic::Ptr(new IceModel_area_glacierized_shelf(this));
  m_ts_diagnostics["dt"]                                   = TSDiagnostic::Ptr(new IceModel_dt(this));
  m_ts_diagnostics["max_diffusivity"]                      = TSDiagnostic::Ptr(new IceModel_max_diffusivity(this));
  m_ts_diagnostics["enthalpy_glacierized"]                 = TSDiagnostic::Ptr(new IceModel_enthalpy_glacierized(this));
  m_ts_diagnostics["enthalpy_nonglacierized"]              = TSDiagnostic::Ptr(new IceModel_enthalpy_nonglacierized(this));
  m_ts_diagnostics["max_hor_vel"]                          = TSDiagnostic::Ptr(new IceModel_max_hor_vel(this));
  m_ts_diagnostics["limnsw"]                               = TSDiagnostic::Ptr(new IceModel_limnsw(this));

  m_ts_diagnostics["surface_ice_flux"]                   = TSDiagnostic::Ptr(new IceModel_surface_flux(this));
  m_ts_diagnostics["surface_ice_flux_cumulative"]        = TSDiagnostic::Ptr(new IceModel_surface_flux_cumulative(this));
  m_ts_diagnostics["grounded_basal_ice_flux"]            = TSDiagnostic::Ptr(new IceModel_grounded_basal_flux(this));
  m_ts_diagnostics["grounded_basal_ice_flux_cumulative"] = TSDiagnostic::Ptr(new IceModel_grounded_basal_flux_cumulative(this));
  m_ts_diagnostics["sub_shelf_ice_flux"]                 = TSDiagnostic::Ptr(new IceModel_sub_shelf_flux(this));
  m_ts_diagnostics["sub_shelf_ice_flux_cumulative"]      = TSDiagnostic::Ptr(new IceModel_sub_shelf_flux_cumulative(this));
  m_ts_diagnostics["nonneg_rule_flux"]                   = TSDiagnostic::Ptr(new IceModel_nonneg_flux(this));
  m_ts_diagnostics["nonneg_rule_flux_cumulative"]        = TSDiagnostic::Ptr(new IceModel_nonneg_flux_cumulative(this));
  m_ts_diagnostics["discharge_flux"]                     = TSDiagnostic::Ptr(new IceModel_discharge_flux(this));
  m_ts_diagnostics["discharge_flux_cumulative"]          = TSDiagnostic::Ptr(new IceModel_discharge_flux_cumulative(this));
  m_ts_diagnostics["H_to_Href_flux"]                     = TSDiagnostic::Ptr(new IceModel_H_to_Href_flux(this));
  m_ts_diagnostics["Href_to_H_flux"]                     = TSDiagnostic::Ptr(new IceModel_Href_to_H_flux(this));
  m_ts_diagnostics["sum_divQ_flux"]                      = TSDiagnostic::Ptr(new IceModel_sum_divQ_flux(this));

  // get diagnostics from submodels
  std::map<std::string, const Component*>::const_iterator j = m_submodels.begin();
  while (j != m_submodels.end()) {
    j->second->get_diagnostics(m_diagnostics, m_ts_diagnostics);
    ++j;
  }
}

void IceModel::list_diagnostics() {

  m_log->message(1, "\n");

  // quantities with dedicated storage
  {
    std::set<std::string> list = m_grid->variables().keys();

    for (unsigned int d = 3; d > 1; --d) {

      m_log->message(1,
                     "======== Available %dD quantities with dedicated storage ========\n",
                     d);

      std::set<std::string>::iterator j;
      for (j = list.begin(); j != list.end(); ++j) {
        const IceModelVec *v = NULL;

        if (m_grid->variables().is_available(*j)) {
          v = m_grid->variables().get(*j);
        }

        if (v != NULL && v->get_ndims() == d) {
          const SpatialVariableMetadata &var = v->metadata();

          std::string
            name                = var.get_name(),
            units               = var.get_string("units"),
            glaciological_units = var.get_string("glaciological_units"),
            long_name           = var.get_string("long_name");

          if (not glaciological_units.empty()) {
            units = glaciological_units;
          }

          m_log->message(1,
                         "   Name: %s [%s]\n"
                         "       - %s\n\n", name.c_str(), units.c_str(), long_name.c_str());
        }
      }
    }

  } // done with quantities with dedicated storage

  // 2D and 3D diagnostics
  for (unsigned int d = 3; d > 1; --d) {

    m_log->message(1,
                   "======== Available %dD diagnostic quantities ========\n",
                   d);

    std::map<std::string, Diagnostic::Ptr>::iterator j = m_diagnostics.begin();
    while (j != m_diagnostics.end()) {
      Diagnostic::Ptr diag = j->second;

      std::string
        name                = j->first,
        units               = diag->get_metadata().get_string("units"),
        glaciological_units = diag->get_metadata().get_string("glaciological_units");

      if (not glaciological_units.empty()) {
        units = glaciological_units;
      }

      if (diag->get_metadata().get_n_spatial_dimensions() == d) {

        m_log->message(1, "   Name: %s [%s]\n", name.c_str(), units.c_str());

        for (int k = 0; k < diag->get_nvars(); ++k) {
          SpatialVariableMetadata var = diag->get_metadata(k);

          std::string long_name = var.get_string("long_name");

          m_log->message(1, "      -  %s\n", long_name.c_str());
        }

        m_log->message(1, "\n");
      }

      ++j;
    }
  }

  // scalar time-series
  m_log->message(1, "======== Available time-series ========\n");

  std::map<std::string, TSDiagnostic::Ptr>::iterator j = m_ts_diagnostics.begin();
  while (j != m_ts_diagnostics.end()) {
    TSDiagnostic::Ptr diag = j->second;

    std::string name = j->first,
      long_name = diag->get_string("long_name"),
      units = diag->get_string("units"),
      glaciological_units = diag->get_string("glaciological_units");

    if (not glaciological_units.empty()) {
      units = glaciological_units;
    }

    m_log->message(1,
                   "   Name: %s [%s]\n"
                   "      -  %s\n\n",
                   name.c_str(), units.c_str(), long_name.c_str());

    ++j;
  }
}


IceModel_hardav::IceModel_hardav(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hardav"));

  // choice to use SSA power; see #285
  const double power = 1.0 / m_config->get_double("stress_balance.ssa.Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("vertical average of ice hardness", "",
            unitstr, unitstr, 0);

  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

//! \brief Computes vertically-averaged ice hardness.
IceModelVec::Ptr IceModel_hardav::compute_impl() {

  const rheology::FlowLaw *flow_law = model->stress_balance()->shallow()->flow_law();
  if (flow_law == NULL) {
    flow_law = model->stress_balance()->modifier()->flow_law();
    if (flow_law == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Can't compute vertically-averaged hardness: no flow law is used.");
    }
  }

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hardav", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2CellType &cell_type = model->cell_type();

  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->ice_thickness();

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(ice_enthalpy);
  list.add(ice_thickness);
  list.add(*result);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *Eij = ice_enthalpy.get_column(i,j);
      const double H = ice_thickness(i,j);
      if (cell_type.icy(i, j)) {
        (*result)(i,j) = rheology::averaged_hardness(*flow_law,
                                                     H, m_grid->kBelowHeight(H),
                                                     &(m_grid->z()[0]), Eij);
      } else { // put negative value below valid range
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


IceModel_rank::IceModel_rank(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "rank"));

  set_attrs("processor rank", "", "1", "", 0);
  m_vars[0].set_time_independent(true);
}

IceModelVec::Ptr IceModel_rank::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "rank", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  IceModelVec::AccessList list;
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    (*result)(p.i(),p.j()) = m_grid->rank();
  }

  return result;
}


IceModel_cts::IceModel_cts(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "cts", m_grid->z()));

  set_attrs("cts = E/E_s(p), so cold-temperate transition surface is at cts = 1", "",
            "", "", 0);
}

IceModelVec::Ptr IceModel_cts::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "cts", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  energy::compute_cts(model->energy_balance_model()->enthalpy(), model->ice_thickness(), *result);

  return result;
}

IceModel_proc_ice_area::IceModel_proc_ice_area(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "proc_ice_area"));

  set_attrs("number of cells containing ice in a processor's domain", "",
            "", "", 0);
  m_vars[0].set_time_independent(true);
}

IceModelVec::Ptr IceModel_proc_ice_area::compute_impl() {

  const IceModelVec2S        &thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &cell_type = model->cell_type();

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "proc_ice_area", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  int ice_filled_cells = 0;

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(thickness);
  for (Points p(*m_grid); p; p.next()) {
    if (cell_type.icy(p.i(), p.j())) {
      ice_filled_cells += 1;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    (*result)(p.i(), p.j()) = ice_filled_cells;
  }

  return result;
}


IceModel_temp::IceModel_temp(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temp", m_grid->z()));

  set_attrs("ice temperature", "land_ice_temperature", "K", "K", 0);
  m_vars[0].set_double("valid_min", 0);
}

IceModelVec::Ptr IceModel_temp::compute_impl() {

  // update vertical levels (in case the grid was extended
  m_vars[0].set_levels(m_grid->z());

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "temp", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3 *enthalpy = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy->get_column(i,j);
      for (unsigned int k=0; k <m_grid->Mz(); ++k) {
        const double depth = (*thickness)(i,j) - m_grid->z(k);
        Tij[k] = EC->temperature(Enthij[k], EC->pressure(depth));
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


IceModel_temp_pa::IceModel_temp_pa(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temp_pa", m_grid->z()));

  set_attrs("pressure-adjusted ice temperature (degrees above pressure-melting point)", "",
            "deg_C", "deg_C", 0);
  m_vars[0].set_double("valid_max", 0);
}

IceModelVec::Ptr IceModel_temp_pa::compute_impl() {
  bool cold_mode = m_config->get_boolean("energy.temperature_based");
  double melting_point_temp = m_config->get_double("constants.fresh_water.melting_point_temperature");

  // update vertical levels (in case the m_grid was extended
  m_vars[0].set_levels(m_grid->z());

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "temp_pa", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3  *enthalpy  = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  double *Tij;
  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Tij = result->get_column(i,j);
      Enthij = enthalpy->get_column(i,j);
      for (unsigned int k=0; k < m_grid->Mz(); ++k) {
        const double depth = (*thickness)(i,j) - m_grid->z(k),
          p = EC->pressure(depth);
        Tij[k] = EC->pressure_adjusted_temperature(Enthij[k], p);

        if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
          // is 273.15
          if (EC->is_temperate_relaxed(Enthij[k],p) && ((*thickness)(i,j) > 0)) {
            Tij[k] = melting_point_temp;
          }
        }

      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

IceModel_temppabase::IceModel_temppabase(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "temppabase"));

  set_attrs("pressure-adjusted ice temperature at the base of ice", "",
            "Celsius", "Celsius", 0);
}

IceModelVec::Ptr IceModel_temppabase::compute_impl() {

  bool cold_mode = m_config->get_boolean("energy.temperature_based");
  double melting_point_temp = m_config->get_double("constants.fresh_water.melting_point_temperature");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "temp_pa_base", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3 *enthalpy = m_grid->variables().get_3d_scalar("enthalpy");

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  const double *Enthij; // columns of these values

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*enthalpy);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points pt(*m_grid); pt; pt.next()) {
      const int i = pt.i(), j = pt.j();

      Enthij = enthalpy->get_column(i,j);

      const double depth = (*thickness)(i,j),
        p = EC->pressure(depth);
      (*result)(i,j) = EC->pressure_adjusted_temperature(Enthij[0], p);

      if (cold_mode) { // if ice is temperate then its pressure-adjusted temp
        // is 273.15
        if (EC->is_temperate_relaxed(Enthij[0],p) && ((*thickness)(i,j) > 0)) {
          (*result)(i,j) = melting_point_temp;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->shift(-melting_point_temp);

  return result;
}

IceModel_enthalpysurf::IceModel_enthalpysurf(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "enthalpysurf"));

  set_attrs("ice enthalpy at 1m below the ice surface", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_enthalpysurf::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "enthalpysurf", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  // compute levels corresponding to 1 m below the ice surface:

  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->ice_thickness();

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) = std::max(ice_thickness(i,j) - 1.0, 0.0);
  }

  ice_enthalpy.getSurfaceValues(*result, *result);  // z=0 slice

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_thickness(i,j) <= 1.0) {
      (*result)(i,j) = m_fill_value;
    }
  }

  return result;
}

IceModel_enthalpybase::IceModel_enthalpybase(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "enthalpybase"));

  set_attrs("ice enthalpy at the base of ice", "",
            "J kg-1", "J kg-1", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_enthalpybase::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "enthalpybase", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->energy_balance_model()->enthalpy().getHorSlice(*result, 0.0);  // z=0 slice

  result->mask_by(model->ice_thickness(), m_fill_value);

  return result;
}


IceModel_tempbase::IceModel_tempbase(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tempbase"));

  set_attrs("ice temperature at the base of ice", "",
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_tempbase::compute_impl() {

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::Ptr enth = IceModel_enthalpybase(model).compute();

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  // result contains basal enthalpy; note that it is allocated by
  // IceModel_enthalpybase::compute().

  const IceModelVec2CellType &cell_type = model->cell_type();

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double depth = (*thickness)(i,j),
        pressure = EC->pressure(depth);
      if (cell_type.icy(i, j)) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}

IceModel_tempsurf::IceModel_tempsurf(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tempsurf"));

  set_attrs("ice temperature at 1m below the ice surface", "",
            "K", "K", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_tempsurf::compute_impl() {

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::Ptr enth = IceModel_enthalpysurf(model).compute();
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(enth);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  // result contains surface enthalpy; note that it is allocated by
  // IceModel_enthalpysurf::compute().

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  double depth = 1.0,
    pressure = EC->pressure(depth);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if ((*thickness)(i,j) > 1) {
        (*result)(i,j) = EC->temperature((*result)(i,j), pressure);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result->metadata(0) = m_vars[0];
  return result;
}


IceModel_liqfrac::IceModel_liqfrac(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "liqfrac", m_grid->z()));

  set_attrs("liquid water fraction in ice (between 0 and 1)", "",
            "1", "1", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("valid_max", 1);
}

IceModelVec::Ptr IceModel_liqfrac::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "liqfrac", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  bool cold_mode = m_config->get_boolean("energy.temperature_based");

  if (cold_mode) {
    result->set(0.0);
  } else {
    energy::compute_liquid_water_fraction(model->energy_balance_model()->enthalpy(),
                                          model->ice_thickness(),
                                          *result);
  }

  return result;
}

IceModel_tempicethk::IceModel_tempicethk(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "tempicethk"));

  set_attrs("temperate ice thickness (total column content)", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_tempicethk::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "tempicethk", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->cell_type();
  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->ice_thickness();

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);
  list.add(ice_enthalpy);
  list.add(ice_thickness);

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        const double *Enth = ice_enthalpy.get_column(i,j);
        double H_temperate = 0.0;
        const double H = ice_thickness(i,j);
        const unsigned int ks = m_grid->kBelowHeight(H);

        for (unsigned int k=0; k<ks; ++k) { // FIXME issue #15
          double pressure = EC->pressure(H - m_grid->z(k));

          if (EC->is_temperate_relaxed(Enth[k], pressure)) {
            H_temperate += m_grid->z(k+1) - m_grid->z(k);
          }
        }

        double pressure = EC->pressure(H - m_grid->z(ks));
        if (EC->is_temperate_relaxed(Enth[ks], pressure)) {
          H_temperate += H - m_grid->z(ks);
        }

        (*result)(i,j) = H_temperate;
      } else {
        // ice-free
        (*result)(i,j) = m_fill_value;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceModel_tempicethk_basal::IceModel_tempicethk_basal(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "tempicethk_basal"));

  set_attrs("thickness of the basal layer of temperate ice", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

/*!
 * Uses linear interpolation to go beyond vertical grid resolution.
 */
IceModelVec::Ptr IceModel_tempicethk_basal::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "tempicethk_basal", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  EnthalpyConverter::Ptr EC = model->ctx()->enthalpy_converter();

  const IceModelVec2CellType &cell_type = model->cell_type();
  const IceModelVec3& ice_enthalpy = model->energy_balance_model()->enthalpy();
  const IceModelVec2S& ice_thickness = model->ice_thickness();

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);
  list.add(ice_thickness);
  list.add(ice_enthalpy);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = ice_thickness(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (cell_type.ice_free(i, j)) {
        (*result)(i,j) = m_fill_value;
        continue;
      }

      const double *Enth = ice_enthalpy.get_column(i,j);

      unsigned int ks = m_grid->kBelowHeight(H);

      unsigned int k = 0;
      double pressure = EC->pressure(H - m_grid->z(k));
      while (k <= ks) {         // FIXME issue #15
        pressure = EC->pressure(H - m_grid->z(k));

        if (EC->is_temperate_relaxed(Enth[k],pressure)) {
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
        (*result)(i,j) = m_grid->z(ks);
        continue;
      }

      double
        pressure_0 = EC->pressure(H - m_grid->z(k-1)),
        dz         = m_grid->z(k) - m_grid->z(k-1),
        slope1     = (Enth[k] - Enth[k-1]) / dz,
        slope2     = (EC->enthalpy_cts(pressure) - EC->enthalpy_cts(pressure_0)) / dz;

      if (slope1 != slope2) {
        (*result)(i,j) = m_grid->z(k-1) +
          (EC->enthalpy_cts(pressure_0) - Enth[k-1]) / (slope1 - slope2);

        // check if the resulting thickness is valid:
        (*result)(i,j) = std::max((*result)(i,j), m_grid->z(k-1));
        (*result)(i,j) = std::min((*result)(i,j), m_grid->z(k));
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Linear interpolation of the thickness of"
                                      " the basal temperate layer failed:\n"
                                      "(i=%d, j=%d, k=%d, ks=%d)\n",
                                      i, j, k, ks);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

IceModel_flux_divergence::IceModel_flux_divergence(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "flux_divergence"));

  set_attrs("flux divergence", "", "m s-1", "m year-1", 0);
}

IceModelVec::Ptr IceModel_flux_divergence::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "flux_divergence", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->flux_divergence());

  return result;
}

IceModel_climatic_mass_balance_cumulative::IceModel_climatic_mass_balance_cumulative(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "climatic_mass_balance_cumulative"));

  set_attrs("cumulative ice-equivalent climatic mass balance", "",
            "kg m-2", "kg m-2", 0);
}

IceModelVec::Ptr IceModel_climatic_mass_balance_cumulative::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "climatic_mass_balance_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->cumulative_fluxes_2d().climatic_mass_balance);

  return result;
}

IceModel_volume_glacierized::IceModel_volume_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);

  m_ts->metadata().set_string("long_name", "volume of the ice in glacierized areas");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_glacierized::update(double a, double b) {

  double value = model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_volume_nonglacierized::IceModel_volume_nonglacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_nonglacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);

  m_ts->metadata().set_string("long_name", "volume of the ice, including seasonal cover");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_nonglacierized::update(double a, double b) {

  double value = model->ice_volume(0.0);

  m_ts->append(value, a, b);
}

IceModel_slvol::IceModel_slvol(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "slvol", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m");
  m_ts->dimension_metadata().set_string("units", m_time_units);

  m_ts->metadata().set_string("long_name", "total sea-level relevant ice IN SEA-LEVEL EQUIVALENT");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_slvol::update(double a, double b) {

  double value = model->sealevel_volume(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_volume_rate_of_change_glacierized::IceModel_volume_rate_of_change_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_rate_of_change_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3 s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->rate_of_change = true;

  m_ts->metadata().set_string("long_name", "rate of change of the ice volume in glacierized areas");
}

void IceModel_volume_rate_of_change_glacierized::update(double a, double b) {

  double value = model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));

  // note that "value" below *should* be the ice volume
  m_ts->append(value, a, b);
}

IceModel_volume_rate_of_change_nonglacierized::IceModel_volume_rate_of_change_nonglacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_rate_of_change_nonglacierized",
                                  m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3 s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->rate_of_change = true;

  m_ts->metadata().set_string("long_name",
                              "rate of change of the ice volume, including seasonal cover");
}

void IceModel_volume_rate_of_change_nonglacierized::update(double a, double b) {

  // note that "value" below *should* be the ice volume
  m_ts->append(model->ice_volume(0.0), a, b);
}


IceModel_area_glacierized::IceModel_area_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "area_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "glacierized area");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_area_glacierized::update(double a, double b) {

  double value = model->ice_area(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_mass_glacierized::IceModel_mass_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "mass_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass of the ice in glacierized areas");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_mass_glacierized::update(double a, double b) {

  double value = model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value * m_grid->ctx()->config()->get_double("constants.ice.density"), a, b);
}

IceModel_mass_nonglacierized::IceModel_mass_nonglacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "mass_nonglacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass of the ice, including seasonal cover");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_mass_nonglacierized::update(double a, double b) {

  double value = model->ice_volume(0.0);

  m_ts->append(value * m_grid->ctx()->config()->get_double("constants.ice.density"), a, b);
}

IceModel_mass_rate_of_change_glacierized::IceModel_mass_rate_of_change_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "mass_rate_of_change_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "rate of change of the mass of ice in glacierized areas");

  m_ts->rate_of_change = true;
}

void IceModel_mass_rate_of_change_glacierized::update(double a, double b) {

  const double
    ice_density = m_grid->ctx()->config()->get_double("constants.ice.density"),
    ice_volume  = model->ice_volume(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(ice_volume * ice_density, a, b);
}

IceModel_mass_rate_of_change_nonglacierized::IceModel_mass_rate_of_change_nonglacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "mass_rate_of_change_nonglacierized",
                                  m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                              "rate of change of the mass of ice, including seasonal cover");

  m_ts->rate_of_change = true;
}

void IceModel_mass_rate_of_change_nonglacierized::update(double a, double b) {

  const double ice_density = m_grid->ctx()->config()->get_double("constants.ice.density");

  m_ts->append(model->ice_volume(0.0) * ice_density, a, b);
}


IceModel_volume_glacierized_temperate::IceModel_volume_glacierized_temperate(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_glacierized_temperate", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of temperate ice in glacierized areas");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_glacierized_temperate::update(double a, double b) {

  double value = model->ice_volume_temperate(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_volume_nonglacierized_temperate::IceModel_volume_nonglacierized_temperate(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_nonglacierized_temperate", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of temperate ice, including seasonal cover");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_nonglacierized_temperate::update(double a, double b) {

  double value = model->ice_volume_temperate(0.0);

  m_ts->append(value, a, b);
}


IceModel_volume_glacierized_cold::IceModel_volume_glacierized_cold(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_glacierized_cold", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of cold ice in glacierized areas");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_glacierized_cold::update(double a, double b) {

  double value = model->ice_volume_cold(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_volume_nonglacierized_cold::IceModel_volume_nonglacierized_cold(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_nonglacierized_cold", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of cold ice, including seasonal cover");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_volume_nonglacierized_cold::update(double a, double b) {

  double value = model->ice_volume_cold(0.0);

  m_ts->append(value, a, b);
}

IceModel_area_glacierized_temperate_base::IceModel_area_glacierized_temperate_base(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "area_glacierized_temperate_base", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "glacierized area where basal ice is temperate");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_area_glacierized_temperate_base::update(double a, double b) {

  double value = model->ice_area_temperate(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_area_glacierized_cold_base::IceModel_area_glacierized_cold_base(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "area_glacierized_cold_base", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "glacierized area where basal ice is cold");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_area_glacierized_cold_base::update(double a, double b) {

  double value = model->ice_area_cold(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_enthalpy_glacierized::IceModel_enthalpy_glacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "enthalpy_glacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "J");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "enthalpy of the ice in glacierized areas");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_enthalpy_glacierized::update(double a, double b) {

  double value = energy::total_ice_enthalpy(m_config->get_double("output.ice_free_thickness_standard"),
                                            model->energy_balance_model()->enthalpy(),
                                            model->ice_thickness(),
                                            model->cell_area());

  m_ts->append(value, a, b);
}

IceModel_enthalpy_nonglacierized::IceModel_enthalpy_nonglacierized(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "enthalpy_nonglacierized", m_time_dimension_name);

  m_ts->metadata().set_string("units", "J");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "enthalpy of the ice, including seasonal cover");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_enthalpy_nonglacierized::update(double a, double b) {

  double value = energy::total_ice_enthalpy(0.0,
                                            model->energy_balance_model()->enthalpy(),
                                            model->ice_thickness(),
                                            model->cell_area());

  m_ts->append(value, a, b);
}

IceModel_area_glacierized_grounded::IceModel_area_glacierized_grounded(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "area_glacierized_grounded", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "area of grounded ice in glacierized areas");
}

void IceModel_area_glacierized_grounded::update(double a, double b) {

  double value = model->ice_area_grounded(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_area_glacierized_shelf::IceModel_area_glacierized_shelf(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "area_glacierized_shelf", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "area of ice shelves in glacierized areas");
}

void IceModel_area_glacierized_shelf::update(double a, double b) {

  double value = model->ice_area_floating(m_config->get_double("output.ice_free_thickness_standard"));

  m_ts->append(value, a, b);
}

IceModel_dt::IceModel_dt(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "dt", m_time_dimension_name);

  m_ts->metadata().set_string("units", "second");
  m_ts->metadata().set_string("glaciological_units", "year");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass continuity time step");
}

void IceModel_dt::update(double a, double b) {

  m_ts->append(model->dt(), a, b);
}

IceModel_max_diffusivity::IceModel_max_diffusivity(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "max_diffusivity", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m2 s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "maximum diffusivity");
}

void IceModel_max_diffusivity::update(double a, double b) {
  double value = model->stress_balance()->max_diffusivity();

  m_ts->append(value, a, b);
}

IceModel_surface_flux::IceModel_surface_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "surface_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total over ice domain of top surface ice mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_surface_flux::update(double a, double b) {

  m_ts->append(model->cumulative_fluxes().surface, a, b);
}

IceModel_surface_flux_cumulative::IceModel_surface_flux_cumulative(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "surface_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total over ice domain of top surface ice mass flux");
}

void IceModel_surface_flux_cumulative::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().surface, a, b);
}

IceModel_grounded_basal_flux::IceModel_grounded_basal_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "grounded_basal_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total over grounded ice domain of basal mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_grounded_basal_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().grounded_basal, a, b);
}

IceModel_grounded_basal_flux_cumulative::IceModel_grounded_basal_flux_cumulative(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "grounded_basal_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total grounded basal mass flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_grounded_basal_flux_cumulative::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().grounded_basal, a, b);
}

IceModel_sub_shelf_flux::IceModel_sub_shelf_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sub_shelf_ice_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "total sub-shelf ice flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_sub_shelf_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().sub_shelf, a, b);
}

IceModel_sub_shelf_flux_cumulative::IceModel_sub_shelf_flux_cumulative(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sub_shelf_ice_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative total sub-shelf ice flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_sub_shelf_flux_cumulative::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().sub_shelf, a, b);
}

IceModel_nonneg_flux::IceModel_nonneg_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "nonneg_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_nonneg_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().nonneg_rule, a, b);
}

IceModel_nonneg_flux_cumulative::IceModel_nonneg_flux_cumulative(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "nonneg_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative 'numerical' ice flux resulting from enforcing the 'thk >= 0' rule");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_nonneg_flux_cumulative::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().nonneg_rule, a, b);
}

IceModel_discharge_flux::IceModel_discharge_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "discharge_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "discharge (calving & icebergs) flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_discharge_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().discharge, a, b);
}

IceModel_discharge_flux_cumulative::IceModel_discharge_flux_cumulative(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "discharge_flux_cumulative", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "cumulative discharge (calving etc.) flux");
  m_ts->metadata().set_string("comment", "positive means ice gain");
}

void IceModel_discharge_flux_cumulative::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().discharge, a, b);
}

IceModel_dHdt::IceModel_dHdt(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "dHdt"));

  set_attrs("ice thickness rate of change", "tendency_of_land_ice_thickness",
            "m s-1", "m year-1", 0);

  double fill_value = units::convert(m_sys, m_fill_value,
                                     "m year-1", "m second-1");

  m_vars[0].set_double("valid_min",  units::convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max",  units::convert(m_sys,  1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_string("cell_methods", "time: mean");

  m_last_ice_thickness.create(m_grid, "last_ice_thickness", WITHOUT_GHOSTS);
  m_last_ice_thickness.set_attrs("internal",
                               "ice thickness at the time of the last report of dHdt",
                               "m", "land_ice_thickness");

  m_last_report_time = GSL_NAN;
}

IceModelVec::Ptr IceModel_dHdt::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "dHdt", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  if (gsl_isnan(m_last_report_time)) {
    result->set(units::convert(m_sys, 2e6, "m year-1", "m second-1"));
  } else {
    const IceModelVec2S& ice_thickness = model->ice_thickness();

    IceModelVec::AccessList list;
    list.add(*result);
    list.add(m_last_ice_thickness);
    list.add(ice_thickness);

    double dt = m_grid->ctx()->time()->current() - m_last_report_time;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (ice_thickness(i, j) - m_last_ice_thickness(i, j)) / dt;
    }
  }

  // Save the ice thickness and the corresponding time:
  this->update_cumulative();

  return result;
}

void IceModel_dHdt::update_cumulative() {
  const IceModelVec2S& ice_thickness = model->ice_thickness();

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(m_last_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    m_last_ice_thickness(i, j) = ice_thickness(i, j);
  }

  m_last_report_time = m_grid->ctx()->time()->current();
}

IceModel_volume_glacierized_grounded::IceModel_volume_glacierized_grounded(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_glacierized_grounded", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of grounded ice in glacierized areas");
}

void IceModel_volume_glacierized_grounded::update(double a, double b) {
  double volume = 0.0;

  const IceModelVec2CellType &cell_type = model->cell_type();

  const IceModelVec2S
    &cell_area     = model->cell_area(),
    &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  const double thickness_threshold = m_config->get_double("output.ice_free_thickness_standard");

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(cell_type);
  list.add(cell_area);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double H = ice_thickness(i, j);

    if (cell_type.grounded(i, j) and H >= thickness_threshold) {
      volume += cell_area(i, j) * H;
    }
  }

  m_ts->append(GlobalSum(m_grid->com, volume), a, b);
}

IceModel_volume_glacierized_shelf::IceModel_volume_glacierized_shelf(IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "volume_glacierized_shelf", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m3");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "volume of ice shelves in glacierized areas");
}

void IceModel_volume_glacierized_shelf::update(double a, double b) {
  double volume = 0.0;

  const IceModelVec2CellType &cell_type = model->cell_type();

  const IceModelVec2S
    &cell_area     = model->cell_area(),
    &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  const double thickness_threshold = m_config->get_double("output.ice_free_thickness_standard");

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(cell_type);
  list.add(cell_area);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double H = ice_thickness(i, j);

    if (cell_type.ocean(i, j) and H >= thickness_threshold) {
      volume += cell_area(i, j) * H;
    }
  }

  m_ts->append(GlobalSum(m_grid->com, volume), a, b);
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
IceModel_max_hor_vel::IceModel_max_hor_vel(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "max_hor_vel", m_time_dimension_name);

  m_ts->metadata().set_string("units", "m second-1");
  m_ts->metadata().set_string("glaciological_units", "m year-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "maximum abs component of horizontal ice velocity"
                                " over grid in last time step during time-series reporting interval");
}

void IceModel_max_hor_vel::update(double a, double b) {

  CFLData cfl = model->stress_balance()->max_timestep_cfl_3d();

  m_ts->append(std::max(cfl.u_max, cfl.v_max), a, b);
}

IceModel_H_to_Href_flux::IceModel_H_to_Href_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "H_to_Href_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass flux from thk to Href");
  m_ts->metadata().set_string("comment", "does not correspond to mass gain or loss");
  m_ts->rate_of_change = true;
}

void IceModel_H_to_Href_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().H_to_Href, a, b);
}

IceModel_Href_to_H_flux::IceModel_Href_to_H_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "Href_to_H_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass flux from Href to thk");
  m_ts->metadata().set_string("comment", "does not correspond to mass gain or loss");
  m_ts->rate_of_change = true;
}

void IceModel_Href_to_H_flux::update(double a, double b) {
  m_ts->append(model->cumulative_fluxes().Href_to_H, a, b);
}

IceModel_sum_divQ_flux::IceModel_sum_divQ_flux(const IceModel *m)
  : TSDiag<IceModel>(m) {
  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "sum_divQ_flux", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "sum(divQ)");
  m_ts->metadata().set_string("comment", "positive means ice gain");
  m_ts->rate_of_change = true;
}

void IceModel_sum_divQ_flux::update(double a, double b) {

  m_ts->append(model->cumulative_fluxes().sum_divQ_SIA +
               model->cumulative_fluxes().sum_divQ_SSA,
               a, b);
}

IceModel_limnsw::IceModel_limnsw(const IceModel *m)
  : TSDiag<IceModel>(m) {

  // set metadata:
  m_ts = new DiagnosticTimeseries(*m_grid, "limnsw", m_time_dimension_name);

  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name", "mass of the ice not displacing sea water");
  m_ts->metadata().set_double("valid_min", 0.0);
}

void IceModel_limnsw::update(double a, double b) {

  const double
    ice_density = m_config->get_double("constants.ice.density"),
    ice_volume  = model->ice_volume_not_displacing_seawater(m_config->get_double("output.ice_free_thickness_standard")),
    ice_mass    = ice_volume * ice_density;

  m_ts->append(ice_mass, a, b);
}

IceModel_nonneg_flux_2D_cumulative::IceModel_nonneg_flux_2D_cumulative(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "nonneg_flux_cumulative"));

  set_attrs("cumulative non-negative rule (thk >= 0) flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_nonneg_flux_2D_cumulative::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->cumulative_fluxes_2d().nonneg);

  return result;
}

IceModel_grounded_basal_flux_2D_cumulative::IceModel_grounded_basal_flux_2D_cumulative(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "grounded_basal_flux_cumulative"));

  set_attrs("cumulative grounded basal flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_grounded_basal_flux_2D_cumulative::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->cumulative_fluxes_2d().basal_grounded);

  return result;
}

IceModel_floating_basal_flux_2D_cumulative::IceModel_floating_basal_flux_2D_cumulative(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "floating_basal_flux_cumulative"));

  set_attrs("cumulative floating basal flux",
            "",                 // no standard name
            "kg m-2", "Gt m-2", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_floating_basal_flux_2D_cumulative::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->cumulative_fluxes_2d().basal_floating);

  return result;
}


IceModel_discharge_flux_2D_cumulative::IceModel_discharge_flux_2D_cumulative(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "discharge_flux_cumulative"));

  set_attrs("cumulative ice discharge (calving) flux",
            "",                 // no standard name
            "kg", "Gt", 0);
  m_vars[0].set_string("comment", "positive means ice gain");
}

IceModelVec::Ptr IceModel_discharge_flux_2D_cumulative::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->cumulative_fluxes_2d().discharge);

  return result;
}

IceModel_discharge_flux_2D::IceModel_discharge_flux_2D(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "discharge_flux"));

  set_attrs("average ice discharge (calving) flux over reporting interval",
            "",                 // no standard name
            "kg second-1", "Gt year-1", 0);
  m_vars[0].set_string("comment", "positive means ice gain");

  double fill_value = units::convert(m_sys, m_fill_value,
                                     "Gt year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_string("cell_methods", "time: mean");

  m_last_cumulative_discharge.create(m_grid, "last_cumulative_discharge", WITHOUT_GHOSTS);
  m_last_cumulative_discharge.set_attrs("internal",
                                        "cumulative discharge at the time of the last report of discharge_flux",
                                        "kg", "");

  m_last_report_time = GSL_NAN;
}

IceModelVec::Ptr IceModel_discharge_flux_2D::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "discharge_flux", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  const IceModelVec2S &cumulative_discharge = model->cumulative_fluxes_2d().discharge;
  const double current_time = m_grid->ctx()->time()->current();

  if (gsl_isnan(m_last_report_time)) {
    const double fill_value = units::convert(m_sys, m_fill_value,
                                             "Gt year-1", "kg second-1");
    result->set(fill_value);
  } else {
    IceModelVec::AccessList list;
    list.add(*result);
    list.add(m_last_cumulative_discharge);
    list.add(cumulative_discharge);

    double dt = current_time - m_last_report_time;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (cumulative_discharge(i, j) - m_last_cumulative_discharge(i, j)) / dt;
    }
  }

  // Save the cumulative discharge and the corresponding time:
  m_last_cumulative_discharge.copy_from(cumulative_discharge);
  m_last_report_time = current_time;

  return result;
}

IceModel_lat_lon_bounds::IceModel_lat_lon_bounds(const IceModel *m,
                                                 const std::string &var_name,
                                                 const std::string &proj_string)
  : Diag<IceModel>(m) {
  assert(var_name == "lat" || var_name == "lon");
  m_var_name = var_name;

  // set metadata:
  std::vector<double> levels(4);
  for (int k = 0; k < 4; ++k) {
    levels[k] = k;
  }

  m_vars.push_back(SpatialVariableMetadata(m_sys, m_var_name + "_bnds", levels));
  m_vars[0].get_z().set_name("nv4");
  m_vars[0].get_z().clear_all_strings();
  m_vars[0].get_z().clear_all_doubles();
  m_vars[0].set_time_independent(true);

  if (m_var_name == "lon") {
    set_attrs("longitude bounds", "", "degree_east", "degree_east", 0);
    m_vars[0].set_double("valid_min", -180);
    m_vars[0].set_double("valid_max", 180);
  } else {
    set_attrs("latitude bounds", "", "degree_north", "degree_north", 0);
    m_vars[0].set_double("valid_min", -90);
    m_vars[0].set_double("valid_max", 90);
  }
  m_vars[0].set_string("coordinates", "");

  m_proj_string = proj_string;

#if (PISM_USE_PROJ4==1)
  // create PROJ.4 objects to check if proj_string is OK.
  Proj lonlat("+proj=latlong +datum=WGS84 +ellps=WGS84");
  Proj pism(m_proj_string);
#endif
  // If PISM_USE_PROJ4 is not 1 we don't need to check validity of m_proj_string: this diagnostic
  // will not be available and so this code will not run.
}

IceModelVec::Ptr IceModel_lat_lon_bounds::compute_impl() {
  std::map<std::string,std::string> attrs;
  std::vector<double> indices(4);

  IceModelVec3Custom::Ptr result(new IceModelVec3Custom);
  result->create(m_grid, m_var_name + "_bnds", "nv4",
                 indices, attrs);
  result->metadata(0) = m_vars[0];

  bool latitude = true;
  if (m_var_name == "lon") {
    latitude = false;
  }

  if (latitude) {
    compute_lat_bounds(m_proj_string, *result);
  } else {
    compute_lon_bounds(m_proj_string, *result);
  }

  return result;
}

IceModel_land_ice_area_fraction::IceModel_land_ice_area_fraction(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, land_ice_area_fraction_name));
  set_attrs("fraction of a grid cell covered by ice (grounded or floating)",
            "",                 // no standard name
            "1", "1", 0);
}

IceModelVec::Ptr IceModel_land_ice_area_fraction::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, land_ice_area_fraction_name, WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const Vars &variables = m_grid->variables();

  const IceModelVec2S
    &thickness         = *variables.get_2d_scalar("land_ice_thickness"),
    &surface_elevation = *variables.get_2d_scalar("surface_altitude"),
    &bed_topography    = *variables.get_2d_scalar("bedrock_altitude");

  const IceModelVec2CellType &cell_type = model->cell_type();

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(surface_elevation);
  list.add(bed_topography);
  list.add(cell_type);
  list.add(*result);

  const bool do_part_grid = m_config->get_boolean("geometry.part_grid.enabled");
  const IceModelVec2S *Href = NULL;
  if (do_part_grid) {
    Href = variables.get_2d_scalar("Href");
    list.add(*Href);
  }

  const bool reduce_frontal_thickness = false;
  const double dx = m_grid->dx();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.icy(i, j)) {
        // an "icy" cell: the area fraction is one
        (*result)(i, j) = 1.0;
      } else if (cell_type.ice_free_ocean(i, j)) {
        // an ice-free ocean cell may be "partially-filled", in which case we need to compute its
        // ice area fraction by dividing Href by the threshold thickness.

        double H_reference = do_part_grid ? (*Href)(i, j) : 0.0;

        if (H_reference > 0.0) {
          const double H_threshold = part_grid_threshold_thickness(cell_type.int_star(i, j),
                                                                   thickness.star(i, j),
                                                                   surface_elevation.star(i, j),
                                                                   bed_topography(i,j),
                                                                   dx,
                                                                   reduce_frontal_thickness);
          // protect from a division by zero
          if (H_threshold > 0.0) {
            (*result)(i, j) = H_reference / H_threshold;
          } else {
            (*result)(i, j) = 1.0;
          }
        } else {
          // H_reference is zero
          (*result)(i, j) = 0.0;
        }
      } else {
        // an ice-free-ground cell: the area fraction is zero
        (*result)(i, j) = 0.0;
      }
    } // end of the loop over grid points
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceModel_grounded_ice_sheet_area_fraction::IceModel_grounded_ice_sheet_area_fraction(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, grounded_ice_sheet_area_fraction_name));
  set_attrs("fraction of a grid cell covered by grounded ice",
            "",                 // no standard name
            "1", "1", 0);
}

IceModelVec::Ptr IceModel_grounded_ice_sheet_area_fraction::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, grounded_ice_sheet_area_fraction_name, WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const double sea_level = model->ocean_model()->sea_level_elevation();

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density");

  const Vars &variables = m_grid->variables();

  const IceModelVec2S
    &ice_thickness  = *variables.get_2d_scalar("land_ice_thickness"),
    &bed_topography = *variables.get_2d_scalar("bedrock_altitude");

  const IceModelVec2CellType &cell_type = model->cell_type();

  compute_grounded_cell_fraction(ice_density, ocean_density, sea_level,
                                 ice_thickness, bed_topography, cell_type,
                                 *result, NULL, NULL);

  // All grounded areas have the grounded cell fraction of one, so now we make sure that ice-free
  // areas get the value of 0 (they are grounded but not covered by a grounded ice sheet).

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if (cell_type.ice_free(i, j)) {
        (*result)(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceModel_floating_ice_sheet_area_fraction::IceModel_floating_ice_sheet_area_fraction(const IceModel *m)
  : Diag<IceModel>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, floating_ice_sheet_area_fraction_name));
  set_attrs("fraction of a grid cell covered by floating ice",
            "",                 // no standard name
            "1", "1", 0);
}

IceModelVec::Ptr IceModel_floating_ice_sheet_area_fraction::compute_impl() {

  IceModel_land_ice_area_fraction land_ice_area_fraction(model);
  IceModelVec::Ptr ice_area_fraction = land_ice_area_fraction.compute();

  IceModel_grounded_ice_sheet_area_fraction grounded_ice_sheet_area_fraction(model);
  IceModelVec::Ptr grounded_area_fraction = grounded_ice_sheet_area_fraction.compute();

  IceModelVec::Ptr result = ice_area_fraction;
  result->metadata() = m_vars[0];

  // Floating area fraction is total area fraction minus grounded area fraction.
  result->add(-1.0, *grounded_area_fraction);

  return result;
}

IceModel_surface_mass_balance_average::IceModel_surface_mass_balance_average(const IceModel *m)
  : Diag<IceModel>(m) {
  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "surface_mass_balance_average"));

  set_attrs("average surface mass flux over reporting interval",
            "",                 // no standard name
            "kg m-2 second-1", "kg m-2 year-1", 0);
  m_vars[0].set_string("comment", "positive means ice gain");

  double fill_value = units::convert(m_sys, m_fill_value,
                                     "kg year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_string("cell_methods", "time: mean");

  m_last_cumulative_SMB.create(m_grid, "last_cumulative_SMB", WITHOUT_GHOSTS);
  m_last_cumulative_SMB.set_attrs("internal",
                                  "cumulative SMB at the time of the last report of surface_mass_balance_average",
                                  "kg m-2", "");

  m_last_report_time = GSL_NAN;
}

IceModelVec::Ptr IceModel_surface_mass_balance_average::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "surface_mass_balance_average", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  const IceModelVec2S &cumulative_SMB = model->cumulative_fluxes_2d().climatic_mass_balance;
  const double current_time = m_grid->ctx()->time()->current();

  if (gsl_isnan(m_last_report_time)) {
    const double fill_value = units::convert(m_sys, m_fill_value,
                                             "kg year-1", "kg second-1");
    result->set(fill_value);
  } else {
    IceModelVec::AccessList list;
    list.add(*result);
    list.add(m_last_cumulative_SMB);
    list.add(cumulative_SMB);

    double dt = current_time - m_last_report_time;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (cumulative_SMB(i, j) - m_last_cumulative_SMB(i, j)) / dt;
    }
  }

  // Save the cumulative SMB and the corresponding time:
  m_last_cumulative_SMB.copy_from(cumulative_SMB);
  m_last_report_time = current_time;

  return result;
}

IceModel_basal_mass_balance_average::IceModel_basal_mass_balance_average(const IceModel *m)
  : Diag<IceModel>(m) {
  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "basal_mass_balance_average"));

  set_attrs("average basal mass flux over reporting interval",
            "",                 // no standard name
            "kg m-2 second-1", "kg m-2 year-1", 0);
  m_vars[0].set_string("comment", "positive means ice gain");

  double fill_value = units::convert(m_sys, m_fill_value,
                                     "kg year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_string("cell_methods", "time: mean");

  m_last_cumulative_BMB.create(m_grid, "last_cumulative_basal_mass_balance", WITHOUT_GHOSTS);
  m_last_cumulative_BMB.set_attrs("internal",
                                  "cumulative basal mass balance at the time of the last report of basal_mass_balance_average",
                                  "kg m-2", "");

  m_last_report_time = GSL_NAN;
}

IceModelVec::Ptr IceModel_basal_mass_balance_average::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "basal_mass_balance_average", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  const IceModelVec2S &cumulative_grounded_BMB = model->cumulative_fluxes_2d().basal_grounded;
  const IceModelVec2S &cumulative_floating_BMB = model->cumulative_fluxes_2d().basal_floating;
  const double current_time = m_grid->ctx()->time()->current();

  if (gsl_isnan(m_last_report_time)) {
    const double fill_value = units::convert(m_sys, m_fill_value,
                                             "kg year-1", "kg second-1");
    result->set(fill_value);
  } else {
    IceModelVec::AccessList list;
    list.add(*result);
    list.add(m_last_cumulative_BMB);
    list.add(cumulative_grounded_BMB);
    list.add(cumulative_floating_BMB);

    double dt = current_time - m_last_report_time;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = (cumulative_grounded_BMB(i, j) +
                         cumulative_floating_BMB(i, j) -
                         m_last_cumulative_BMB(i, j)) / dt;
    }
  }

  // Save the cumulative BMB and the corresponding time:
  cumulative_grounded_BMB.add(1.0, cumulative_floating_BMB, m_last_cumulative_BMB);
  m_last_report_time = current_time;

  return result;
}

IceModel_height_above_flotation::IceModel_height_above_flotation(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys,
                                           "height_above_flotation"));

  set_attrs("the height above flotation", "",
            "m", "m", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_height_above_flotation::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "height_above_flotation", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->cell_type();

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    ocean_density = m_config->get_double("constants.sea_water.density"),
    sea_level = model->ocean_model()->sea_level_elevation();

  const Vars &variables = m_grid->variables();

  const IceModelVec2S
    &ice_thickness  = *variables.get_2d_scalar("land_ice_thickness"),
    &bed_topography = *variables.get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);
  list.add(ice_thickness);
  list.add(bed_topography);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double thk = ice_thickness(i,j);
      double bed = bed_topography(i,j);

      // if we have no ice, go on to the next grid point (this cell will be
      // marked as "missing" later)
      if (cell_type.ice_free(i, j)) {
        (*result)(i,j) = m_fill_value;
        continue;
      }
      (*result)(i,j) = ((ice_density / ocean_density) * thk) + (bed - sea_level);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}

IceModel_ice_mass::IceModel_ice_mass(const IceModel *m)
  : Diag<IceModel>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "ice_mass"));

  set_attrs("mass per cell",
            "",                 // no standard name
            "kg", "kg", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

IceModelVec::Ptr IceModel_ice_mass::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "ice_mass", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2CellType &cell_type = model->cell_type();

  const double
    ice_density = m_config->get_double("constants.ice.density");

  const Vars &variables = m_grid->variables();

  const IceModelVec2S
    &ice_thickness = *variables.get_2d_scalar("land_ice_thickness"),
    &cell_area     = *variables.get_2d_scalar("cell_area");

  IceModelVec::AccessList list;
  list.add(cell_type);
  list.add(*result);
  list.add(ice_thickness);
  list.add(cell_area);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // count all ice, including cells which have so little they
      // are considered "ice-free"
      if (ice_thickness(i, j) > 0.0) {
        (*result)(i,j) = ice_density * ice_thickness(i, j) * cell_area(i, j);
      } else {
        (*result)(i,j) = m_fill_value;
      }
    } // end of loop over grid points

  } catch (...) {
    loop.failed();
  }
  loop.check();

  // Add the mass of ice in Href:
  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    const IceModelVec2S &Href = *variables.get_2d_scalar("Href");
    list.add(Href);
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (ice_thickness(i, j) <= 0.0 and Href(i, j) > 0.0) {
        (*result)(i, j) = ice_density * Href(i, j) * cell_area(i,j);
      }
    }
  }

  return result;
}

IceModel_topg_sl_adjusted::IceModel_topg_sl_adjusted(const IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars.push_back(SpatialVariableMetadata(m_sys, "topg_sl_adjusted"));

  set_attrs("sea-level adjusted bed topography (zero is at sea level)", "",
            "meters", "meters", 0);
}

IceModelVec::Ptr IceModel_topg_sl_adjusted::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "topg_sl_adjusted", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->bed_model()->bed_elevation());
  // result = topg - sea_level
  result->shift(-model->ocean_model()->sea_level_elevation());

  return result;
}

IceModel_hardness::IceModel_hardness(const IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hardness", m_grid->z()));

  const double power = 1.0 / m_config->get_double("stress_balance.sia.Glen_exponent");
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);

  set_attrs("ice hardness computed using the SIA flow law", "",
            unitstr, unitstr, 0);
}

IceModelVec::Ptr IceModel_hardness::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "hardness", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const IceModelVec3  &ice_enthalpy  = model->energy_balance_model()->enthalpy();
  const IceModelVec2S &ice_thickness = model->ice_thickness();

  const rheology::FlowLaw *flow_law = model->stress_balance()->modifier()->flow_law();

  IceModelVec::AccessList list;
  list.add(ice_enthalpy);
  list.add(ice_thickness);
  list.add(*result);

  const unsigned int Mz = m_grid->Mz();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      const double *E = ice_enthalpy.get_column(i, j);
      const double H = ice_thickness(i, j);

      double *hardness = result->get_column(i, j);

      for (unsigned int k = 0; k < Mz; ++k) {
        const double depth = H - m_grid->z(k);

        // EC->pressure() handles negative depths correctly
        const double pressure = EC->pressure(depth);

        hardness[k] = flow_law->hardness(E[k], pressure);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

IceModel_viscosity::IceModel_viscosity(IceModel *m)
  : Diag<IceModel>(m) {

  /* set metadata: */
  m_vars.push_back(SpatialVariableMetadata(m_sys, "effective_viscosity",
                                           m_grid->z()));

  set_attrs("effective viscosity of ice", "",
            "Pascal second", "kPascal second", 0);
  m_vars[0].set_double("valid_min", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}

static inline double square(double x) {
  return x * x;
}

IceModelVec::Ptr IceModel_viscosity::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "effective_viscosity", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  IceModelVec3 W;
  W.create(m_grid, "wvel", WITH_GHOSTS);

  using mask::ice_free;

  EnthalpyConverter::Ptr EC = m_grid->ctx()->enthalpy_converter();

  const rheology::FlowLaw *flow_law = model->stress_balance()->modifier()->flow_law();

  const IceModelVec2S &ice_thickness = model->ice_thickness();

  const IceModelVec3
    &ice_enthalpy     = model->energy_balance_model()->enthalpy(),
    &U                = model->stress_balance()->velocity_u(),
    &V                = model->stress_balance()->velocity_v(),
    &W_without_ghosts = model->stress_balance()->velocity_w();

  W_without_ghosts.update_ghosts(W);

  const unsigned int Mz = m_grid->Mz();
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();
  const std::vector<double> &z = m_grid->z();

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(U);
  list.add(V);
  list.add(W);
  list.add(ice_enthalpy);
  list.add(ice_thickness);
  list.add(mask);
  list.add(*result);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double *E = ice_enthalpy.get_column(i, j);
      const double H = ice_thickness(i, j);

      const double
        *u   = U.get_column(i, j),
        *u_n = U.get_column(i, j + 1),
        *u_e = U.get_column(i + 1, j),
        *u_s = U.get_column(i, j - 1),
        *u_w = U.get_column(i - 1, j);

      const double
        *v   = V.get_column(i, j),
        *v_n = V.get_column(i, j + 1),
        *v_e = V.get_column(i + 1, j),
        *v_s = V.get_column(i, j - 1),
        *v_w = V.get_column(i - 1, j);

      const double
        *w   = W.get_column(i, j),
        *w_n = W.get_column(i, j + 1),
        *w_e = W.get_column(i + 1, j),
        *w_s = W.get_column(i, j - 1),
        *w_w = W.get_column(i - 1, j);

      StarStencil<int> m = mask.int_star(i, j);
      const unsigned int
        east  = ice_free(m.e) ? 0 : 1,
        west  = ice_free(m.w) ? 0 : 1,
        south = ice_free(m.s) ? 0 : 1,
        north = ice_free(m.n) ? 0 : 1;

      double *viscosity = result->get_column(i, j);

      if (ice_free(m.ij)) {
        result->set_column(i, j, m_fill_value);
        continue;
      }

      for (unsigned int k = 0; k < Mz; ++k) {
        const double depth = H - z[k];

        if (depth < 0.0) {
          viscosity[k] = m_fill_value;
          continue;
        }

        // EC->pressure() handles negative depths correctly
        const double pressure = EC->pressure(depth);

        const double hardness = flow_law->hardness(E[k], pressure);

        double u_x = 0.0, v_x = 0.0, w_x = 0.0;
        if (west + east > 0) {
          const double D = 1.0 / (dx * (west + east));
          u_x = D * (west * (u[k] - u_w[k]) + east * (u_e[k] - u[k]));
          v_x = D * (west * (v[k] - v_w[k]) + east * (v_e[k] - v[k]));
          w_x = D * (west * (w[k] - w_w[k]) + east * (w_e[k] - w[k]));
        }

        double u_y = 0.0, v_y = 0.0, w_y = 0.0;
        if (south + north > 0) {
          const double D = 1.0 / (dy * (south + north));
          u_y = D * (south * (u[k] - u_s[k]) + north * (u_n[k] - u[k]));
          v_y = D * (south * (v[k] - v_s[k]) + north * (v_n[k] - v[k]));
          w_y = D * (south * (w[k] - w_s[k]) + north * (w_n[k] - w[k]));
        }

        double
          u_z = 0.0,
          v_z = 0.0,
          w_z = 0.0;

        if (k == 0) {
          const double dz = z[1] - z[0];
          u_z = (u[1] - u[0]) / dz;
          v_z = (v[1] - v[0]) / dz;
          w_z = (w[1] - w[0]) / dz;
        } else if (k == Mz - 1) {
          const double dz = z[Mz - 1] - z[Mz - 2];
          u_z = (u[Mz - 1] - u[Mz - 2]) / dz;
          v_z = (v[Mz - 1] - v[Mz - 2]) / dz;
          w_z = (w[Mz - 1] - w[Mz - 2]) / dz;
        } else {
          const double
            dz_p = z[k + 1] - z[k],
            dz_m = z[k] - z[k - 1];
          u_z = 0.5 * ((u[k + 1] - u[k]) / dz_p + (u[k] - u[k - 1]) / dz_m);
          v_z = 0.5 * ((v[k + 1] - v[k]) / dz_p + (v[k] - v[k - 1]) / dz_m);
          w_z = 0.5 * ((w[k + 1] - w[k]) / dz_p + (w[k] - w[k - 1]) / dz_m);
        }

        // These should be "epsilon dot", but that's just too long.
        const double
          eps_xx = u_x,
          eps_yy = v_y,
          eps_zz = w_z,
          eps_xy = 0.5 * (u_y + v_x),
          eps_xz = 0.5 * (u_z + w_x),
          eps_yz = 0.5 * (v_z + w_y);

        // The second invariant of the 3D strain rate tensor; see equation 4.8 in [@ref
        // GreveBlatter2009]. Unlike secondInvariant_2D(), this code does not make assumptions about
        // the input velocity field: we do not ignore w_x and w_y and do not assume that u_z and v_z
        // are zero.
        const double
          gamma = (square(eps_xx) + square(eps_yy) + square(eps_zz) +
                   2.0 * (square(eps_xy) + square(eps_xz) + square(eps_yz)));

        double nu = 0.0;
        // Note: in PISM gamma has an extra factor of 1/2; compare to
        flow_law->effective_viscosity(hardness, 0.5 * gamma, &nu, NULL);

        viscosity[k] = nu;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}



} // end of namespace pism
