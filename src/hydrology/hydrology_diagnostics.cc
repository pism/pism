// Copyright (C) 2012-2017 PISM Authors
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
#include "pism/util/Vars.hh"

namespace pism {
namespace hydrology {

BasalWaterThickness::BasalWaterThickness(const Routing *m)
  : Diag<Routing>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwat")};
  set_attrs("thickness of transportable water in subglacial layer", "", "m", "m", 0);
}

IceModelVec::Ptr BasalWaterThickness::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwat", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  result->copy_from(model->subglacial_water_thickness());
  return result;
}

BasalWaterPressure::BasalWaterPressure(const Routing *m)
  : Diag<Routing>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwp")};
  set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
}


IceModelVec::Ptr BasalWaterPressure::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwp", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];
  result->copy_from(model->subglacial_water_pressure());
  return result;
}


RelativeBasalWaterPressure::RelativeBasalWaterPressure(const Routing *m)
  : Diag<Routing>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "bwprel")};
  set_attrs("pressure of transportable water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  m_vars[0].set_double("_FillValue", m_fill_value);
}


IceModelVec::Ptr RelativeBasalWaterPressure::compute_impl() const {
  double fill_value = m_fill_value;

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwprel", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec2S
    &P  = model->subglacial_water_pressure(),
    &Po = model->overburden_pressure();

  IceModelVec::AccessList list{result.get(), &Po, &P};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Po(i,j) > 0.0) {
      (*result)(i,j) = P(i, j) / Po(i,j);
    } else {
      (*result)(i,j) = fill_value;
    }
  }

  return result;
}


EffectiveBasalWaterPressure::EffectiveBasalWaterPressure(const Routing *m)
  : Diag<Routing>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "effbwp")};
  set_attrs("effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


IceModelVec::Ptr EffectiveBasalWaterPressure::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "effbwp", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2S
    &P  = model->subglacial_water_pressure(),
    &Po = model->overburden_pressure();

  IceModelVec::AccessList list{&Po, &P, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i, j) = Po(i, j) - P(i, j);
  }

  return result;
}

WallMelt::WallMelt(const Routing *m)
  : Diag<Routing>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "wallmelt")};
  set_attrs("wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr WallMelt::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "wallmelt", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  const IceModelVec2S &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  wall_melt(*model, bed_elevation, *result);
  return result;
}

} // end of namespace hydrology
} // end of namespace pism
