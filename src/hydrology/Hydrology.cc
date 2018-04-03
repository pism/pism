// Copyright (C) 2012-2018 PISM Authors
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

#include "Hydrology.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

namespace diagnostics {

/*! @brief Report total input rate of subglacial water (basal melt rate plus input from
  the surface).
 */
class TotalInputRate : public DiagAverageRate<Hydrology>
{
public:
  TotalInputRate(const Hydrology *m)
    : DiagAverageRate<Hydrology>(m, "subglacial_water_input_rate", RATE) {

    m_vars = {SpatialVariableMetadata(m_sys, "subglacial_water_input_rate")};
    m_accumulator.metadata().set_string("units", "m");

    set_attrs("total input rate of subglacial water "
              "(basal melt rate plus input from the surface)", "",
              "m second-1", "m year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
    m_vars[0].set_string("comment", "positive flux corresponds to water gain");
  }

protected:
  const IceModelVec2S& model_input() {
    return model->total_input_rate();
  }
};

} // end of namespace diagnostics

Inputs::Inputs() {
  surface_input_rate = NULL;
  basal_melt_rate    = NULL;
  cell_type          = NULL;
  ice_thickness      = NULL;
  bed_elevation      = NULL;
  ice_sliding_speed  = NULL;
  no_model_mask      = NULL;
}

Hydrology::Hydrology(IceGrid::ConstPtr g)
  : Component(g) {
  m_hold_bmelt = false;

  m_input_rate.create(m_grid, "total_input", WITHOUT_GHOSTS);
  m_input_rate.set_attrs("internal",
                         "hydrology model workspace for total input rate into subglacial water layer",
                         "m s-1", "");

  // *all* Hydrology classes have layer of water stored in till as a state variable
  m_Wtill.create(m_grid, "tillwat", WITHOUT_GHOSTS);
  m_Wtill.set_attrs("model_state",
                    "effective thickness of subglacial water stored in till",
                    "m", "");
  m_Wtill.metadata().set_double("valid_min", 0.0);

  m_Pover.create(m_grid, "overburden_pressure", WITHOUT_GHOSTS);
  m_Pover.set_attrs("internal", "overburden pressure", "Pa", "");
  m_Pover.metadata().set_double("valid_min", 0.0);

  // needs ghosts in Routing and Distributed
  m_W.create(m_grid, "bwat", WITH_GHOSTS, 1);
  m_W.set_attrs("diagnostic",
                "thickness of transportable subglacial water layer",
                "m", "");
  m_W.metadata().set_double("valid_min", 0.0);
}


Hydrology::~Hydrology() {
  // empty
}

void Hydrology::restart(const PIO &input_file, int record) {
  initialization_message();
  this->restart_impl(input_file, record);
}

void Hydrology::bootstrap(const PIO &input_file,
                          const IceModelVec2S &ice_thickness) {
  initialization_message();
  this->bootstrap_impl(input_file, ice_thickness);
}

void Hydrology::initialize(const IceModelVec2S &W_till,
                           const IceModelVec2S &W,
                           const IceModelVec2S &P) {
  initialization_message();
  this->initialize_impl(W_till, W, P);
}

void Hydrology::restart_impl(const PIO &input_file, int record) {
  m_Wtill.read(input_file, record);

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_Wtill);
}

void Hydrology::bootstrap_impl(const PIO &input_file,
                               const IceModelVec2S &ice_thickness) {
  (void) ice_thickness;

  double tillwat_default = m_config->get_double("bootstrapping.defaults.tillwat");
  m_Wtill.regrid(input_file, OPTIONAL, tillwat_default);

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_Wtill);
}

void Hydrology::initialize_impl(const IceModelVec2S &W_till,
                                const IceModelVec2S &W,
                                const IceModelVec2S &P) {
  (void) W;
  (void) P;
  m_Wtill.copy_from(W_till);
}

void Hydrology::update(double t, double dt, const Inputs& inputs) {
  this->update_impl(t, dt, inputs);
}

std::map<std::string, Diagnostic::Ptr> Hydrology::diagnostics_impl() const {
  using namespace diagnostics;
  return {{"tillwat", Diagnostic::wrap(m_Wtill)},
      {"subglacial_water_input_rate", Diagnostic::Ptr(new TotalInputRate(this))}};
}

void Hydrology::define_model_state_impl(const PIO &output) const {
  m_Wtill.define(output);
}

void Hydrology::write_model_state_impl(const PIO &output) const {
  m_Wtill.write(output);
}

//! Update the overburden pressure from ice thickness.
/*!
  Uses the standard hydrostatic (shallow) approximation of overburden pressure,
  \f[ P_0 = \rho_i g H \f]
*/
void Hydrology::compute_overburden_pressure(const IceModelVec2S &ice_thickness,
                                            IceModelVec2S &result) const {
  // FIXME issue #15

  const double
    ice_density      = m_config->get_double("constants.ice.density"),
    standard_gravity = m_config->get_double("constants.standard_gravity");

  IceModelVec::AccessList list{&ice_thickness, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = ice_density * standard_gravity * ice_thickness(i, j);
  }
}

const IceModelVec2S& Hydrology::overburden_pressure() const {
  return m_Pover;
}

//! Return the effective thickness of the water stored in till.
const IceModelVec2S& Hydrology::till_water_thickness() const {
  return m_Wtill;
}

//! Return the effective thickness of the transportable basal water layer.
const IceModelVec2S& Hydrology::subglacial_water_thickness() const {
  return m_W;
}

const IceModelVec2S& Hydrology::total_input_rate() const {
  return m_input_rate;
}

/*!
  Checks @f$ 0 \le W \le W^{max} @f$.
*/
void check_bounds(const IceModelVec2S& W, double W_max) {

  std::string name = W.metadata().get_string("long_name");

  IceGrid::ConstPtr grid = W.grid();

  IceModelVec::AccessList list(W);
  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (W(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: negative %s of %.6f m at (i, j) = (%d, %d)",
                                      name.c_str(), W(i, j), i, j);
      }

      if (W(i,j) > W_max) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: %s of %.6f m exceeds the maximum value %.6f at (i, j) = (%d, %d)",
                                      name.c_str(), W(i, j), W_max, i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Compute the total water input rate into the basal hydrology layer in the ice-covered
//! region.
/*!
  This method ignores the input rate in the ice-free region.

  @param[in] mask cell type mask
  @param[in] basal melt rate (ice thickness per time)
  @param[in] surface_input_rate surface input rate (water thickness per time); set to NULL to ignore
  @param[out] result resulting input rate (water thickness per time)
*/
void Hydrology::compute_input_rate(const IceModelVec2CellType &mask,
                                   const IceModelVec2S &basal_melt_rate,
                                   const IceModelVec2S *surface_input_rate,
                                   IceModelVec2S &result) {

  IceModelVec::AccessList list{&basal_melt_rate, &mask, &result};

  if (surface_input_rate) {
    list.add(*surface_input_rate);
  }

  const double
    ice_density   = m_config->get_double("constants.ice.density"),
    water_density = m_config->get_double("constants.fresh_water.density"),
    C             = ice_density / water_density;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {

      double surface_input = surface_input_rate ? (*surface_input_rate)(i, j) : 0.0;

      result(i,j) = C * basal_melt_rate(i, j) + surface_input;
    } else {
      result(i,j) = 0.0;
    }
  }
}


} // end of namespace hydrology
} // end of namespace pism
