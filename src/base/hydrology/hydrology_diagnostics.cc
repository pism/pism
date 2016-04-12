// Copyright (C) 2012-2016 PISM Authors
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


#include "hydrology_diagnostics.hh"

namespace pism {
namespace hydrology {

Hydrology_bwat::Hydrology_bwat(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "bwat"));
  set_attrs("thickness of transportable water in subglacial layer", "", "m", "m", 0);
}

IceModelVec::Ptr Hydrology_bwat::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "bwat", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;
  model->subglacial_water_thickness(*result);
  return result;
}

Hydrology_bwp::Hydrology_bwp(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "bwp"));
  set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
}


IceModelVec::Ptr Hydrology_bwp::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "bwp", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;
  model->subglacial_water_pressure(*result);
  return result;
}


Hydrology_bwprel::Hydrology_bwprel(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "bwprel"));
  set_attrs("pressure of transportable water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  m_vars[0].set_double("_FillValue", m_config->get_double("fill_value"));
}


IceModelVec::Ptr Hydrology_bwprel::compute_impl() {
  double fill_value = m_config->get_double("fill_value");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "bwprel", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  IceModelVec2S Po;
  Po.create(m_grid, "Po_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(*result);
  model->overburden_pressure(Po);

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(Po);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Po(i,j) > 0.0) {
      (*result)(i,j) /= Po(i,j);
    } else {
      (*result)(i,j) = fill_value;
    }
  }

  return result;
}


Hydrology_effbwp::Hydrology_effbwp(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "effbwp"));
  set_attrs("effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


IceModelVec::Ptr Hydrology_effbwp::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "effbwp", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  IceModelVec2S P;
  P.create(m_grid, "P_temporary", WITHOUT_GHOSTS);

  model->subglacial_water_pressure(P);
  model->overburden_pressure(*result);
  result->add(-1.0, P);  // result <-- result + (-1.0) P = Po - P

  return result;
}


Hydrology_hydrobmelt::Hydrology_hydrobmelt(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hydrobmelt"));
  set_attrs("the version of bmelt seen by the hydrology model",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_hydrobmelt::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hydrobmelt", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->write_in_glaciological_units = true;

  // the value reported diagnostically is merely the last value filled
  result->copy_from(model->m_bmelt_local);

  return result;
}


Hydrology_hydroinput::Hydrology_hydroinput(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "hydroinput"));
  set_attrs("total water input into subglacial hydrology layer",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_hydroinput::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hydroinput", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  // the value reported diagnostically is merely the last value filled
  result->copy_from(model->m_total_input);

  return result;
}


Hydrology_wallmelt::Hydrology_wallmelt(Hydrology *m)
  : Diag<Hydrology>(m) {
  m_vars.push_back(SpatialVariableMetadata(m_sys, "wallmelt"));
  set_attrs("wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Hydrology_wallmelt::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "wallmelt", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  model->wall_melt(*result);

  return result;
}


MCHydrology_ice_free_land_loss_cumulative::MCHydrology_ice_free_land_loss_cumulative(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_ice_free_land_loss_cumulative", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "cumulative liquid water loss from subglacial hydrology into cells with mask as ice free land");
}

void MCHydrology_ice_free_land_loss_cumulative::update(double a, double b) {
  m_ts->append(model->m_ice_free_land_loss_cumulative, a, b);
}

MCHydrology_ice_free_land_loss::MCHydrology_ice_free_land_loss(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_ice_free_land_loss", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "rate of liquid water loss from subglacial hydrology into cells with mask as ice free land");
  m_ts->rate_of_change = true;
}

void MCHydrology_ice_free_land_loss::update(double a, double b) {
  m_ts->append(model->m_ice_free_land_loss_cumulative, a, b);
}

MCHydrology_ocean_loss_cumulative::MCHydrology_ocean_loss_cumulative(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_ocean_loss_cumulative", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "cumulative liquid water loss from subglacial hydrology into cells with mask as ocean");
}

void MCHydrology_ocean_loss_cumulative::update(double a, double b) {
  m_ts->append(model->m_ocean_loss_cumulative, a, b);
}

MCHydrology_ocean_loss::MCHydrology_ocean_loss(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_ocean_loss", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "rate of liquid water loss from subglacial hydrology into cells with mask as ocean");
  m_ts->rate_of_change = true;
}

void MCHydrology_ocean_loss::update(double a, double b) {
  m_ts->append(model->m_ocean_loss_cumulative, a, b);
}

MCHydrology_negative_thickness_gain_cumulative::MCHydrology_negative_thickness_gain_cumulative(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_negative_thickness_gain_cumulative", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "cumulative non-conserving liquid water gain from subglacial hydrology transportable water thickness coming out negative during time step, and being projected up to zero");
}

void MCHydrology_negative_thickness_gain_cumulative::update(double a, double b) {
  m_ts->append(model->m_negative_thickness_gain_cumulative, a, b);
}

MCHydrology_negative_thickness_gain::MCHydrology_negative_thickness_gain(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_negative_thickness_gain", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "rate of non-conserving liquid water gain from subglacial hydrology transportable water thickness coming out negative during time step, and being projected up to zero");
  m_ts->rate_of_change = true;
}

void MCHydrology_negative_thickness_gain::update(double a, double b) {
  m_ts->append(model->m_negative_thickness_gain_cumulative, a, b);
}

MCHydrology_null_strip_loss_cumulative::MCHydrology_null_strip_loss_cumulative(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_null_strip_loss_cumulative", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "cumulative liquid water loss from subglacial hydrology into cells inside the null strip");
}

void MCHydrology_null_strip_loss_cumulative::update(double a, double b) {
  m_ts->append(model->m_null_strip_loss_cumulative, a, b);
}

MCHydrology_null_strip_loss::MCHydrology_null_strip_loss(Routing *m)
      : TSDiag<Routing>(m) {
  m_ts = new DiagnosticTimeseries(*m_grid, "hydro_null_strip_loss", m_time_dimension_name);
  m_ts->metadata().set_string("units", "kg s-1");
  m_ts->dimension_metadata().set_string("units", m_time_units);
  m_ts->metadata().set_string("long_name",
                                "rate of liquid water loss from subglacial hydrology into cells inside the null strip");
  m_ts->rate_of_change = true;
}

void MCHydrology_null_strip_loss::update(double a, double b) {
  m_ts->append(model->m_null_strip_loss_cumulative, a, b);
}

} // end of namespace hydrology
} // end of namespace pism
