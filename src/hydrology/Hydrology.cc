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
}

Hydrology::Hydrology(IceGrid::ConstPtr g)
  : Component(g) {
  m_hold_bmelt = false;

  m_input_rate.create(m_grid, "total_input", WITHOUT_GHOSTS);
  m_input_rate.set_attrs("internal",
                         "hydrology model workspace for total input rate into subglacial water layer",
                         "m s-1", "");

  m_bmelt_local.create(m_grid, "basal_melt_rate_grounded", WITHOUT_GHOSTS);
  m_bmelt_local.set_attrs("internal",
                          "hydrology model workspace for basal_melt_rate_grounded",
                          "m s-1", "");

  // *all* Hydrology classes have layer of water stored in till as a state variable
  m_Wtil.create(m_grid, "tillwat", WITHOUT_GHOSTS);
  m_Wtil.set_attrs("model_state",
                   "effective thickness of subglacial water stored in till",
                   "m", "");
  m_Wtil.metadata().set_double("valid_min", 0.0);

  m_Pover.create(m_grid, "overburden_pressure", WITHOUT_GHOSTS);
  m_Pover.set_attrs("internal", "overburden pressure", "Pa", "");
  m_Pover.metadata().set_double("valid_min", 0.0);
}


Hydrology::~Hydrology() {
  // empty
}


void Hydrology::init() {

  options::String bmelt_file("-hydrology_bmelt_file",
                             "Read time-independent values for basal_melt_rate_grounded from a file;"
                             " replaces basal_melt_rate_grounded computed through conservation of energy");
  // itb = input_to_bed
  options::String itb_file("-hydrology_input_to_bed_file",
                           "A time- and space-dependent file with amount of water"
                           " (depth per time) which should be added to the amount of water"
                           " at the ice sheet bed at the given location at the given time;"
                           " adds to basal_melt_rate_grounded");

  options::Real itb_period_years("-hydrology_input_to_bed_period",
                                 "The period (i.e. duration before repeat), in years,"
                                 " of -hydrology_input_to_bed_file data", 0.0);

  options::Real itb_reference_year("-hydrology_input_to_bed_reference_year",
                                   "The reference year for periodizing the"
                                   " -hydrology_input_to_bed_file data", 0.0);

  // the following are IceModelVec pointers into IceModel generally and are read by code in the
  // update() method at the current Hydrology time

  if (bmelt_file.is_set()) {
    m_log->message(2,
               "  option -hydrology_bmelt_file seen; reading basal_melt_rate_grounded from '%s'.\n", bmelt_file->c_str());
    m_bmelt_local.regrid(bmelt_file, CRITICAL);
    m_hold_bmelt = true;
  }

  InputOptions opts = process_input_options(m_grid->com);

  double tillwat_default = m_config->get_double("bootstrapping.defaults.tillwat");

  switch (opts.type) {
  case INIT_RESTART:
    m_Wtil.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_Wtil.regrid(opts.filename, OPTIONAL, tillwat_default);
    break;
  case INIT_OTHER:
  default:
    m_Wtil.set(tillwat_default);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_Wtil);
}

void Hydrology::update(double t, double dt, const Inputs& inputs) {
  this->update_impl(t, dt, inputs);
}

std::map<std::string, Diagnostic::Ptr> Hydrology::diagnostics_impl() const {
  return {{"tillwat", Diagnostic::wrap(m_Wtil)}};
}

void Hydrology::define_model_state_impl(const PIO &output) const {
  m_Wtil.define(output);
}

void Hydrology::write_model_state_impl(const PIO &output) const {
  m_Wtil.write(output);
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

  IceModelVec::AccessList list{&ice_thickness, &m_Pover};

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
  return m_Wtil;
}

/*!
Checks \f$0 \le W_{til} \le W_{til}^{max} =\f$hydrology_tillwat_max.
 */
void Hydrology::check_Wtil_bounds() {
  double tillwat_max = m_config->get_double("hydrology.tillwat_max");

  IceModelVec::AccessList list(m_Wtil);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_Wtil(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: negative till water effective layer thickness"
                                      " Wtil(i,j) = %.6f m\n"
                                      "at (i,j)=(%d,%d)", m_Wtil(i,j), i, j);
      }

      if (m_Wtil(i,j) > tillwat_max) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "Hydrology: till water effective layer thickness"
                                      " Wtil(i,j) = %.6f m exceeds\n"
                                      "hydrology_tillwat_max = %.6f at (i,j)=(%d,%d)",
                                      m_Wtil(i,j), tillwat_max, i, j);
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
