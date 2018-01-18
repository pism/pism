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
  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t, my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(t - m_t) < 1e-12) && (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  // update Component times: t = current time, t+dt = target time
  m_t  = t;
  m_dt = dt;

  this->update_impl(t, dt, inputs);
}

std::map<std::string, Diagnostic::Ptr> Hydrology::diagnostics_impl() const {
  return {{"tillwat", Diagnostic::wrap(m_Wtill)}};
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

/*!
  Checks \f$0 \le W_{till} \le W_{till}^{max} =\f$hydrology_tillwat_max.
*/
void Hydrology::check_Wtill_bounds() {
  double tillwat_max = m_config->get_double("hydrology.tillwat_max");

  IceModelVec::AccessList list(m_Wtill);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_Wtill(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: negative till water effective layer thickness"
                                      " Wtill(i,j) = %.6f m\n"
                                      "at (i,j)=(%d,%d)", m_Wtill(i,j), i, j);
      }

      if (m_Wtill(i,j) > tillwat_max) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: till water effective layer thickness"
                                      " Wtill(i,j) = %.6f m exceeds\n"
                                      "hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)",
                                      m_Wtill(i,j), tillwat_max, i, j);
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
