// Copyright (C) 2018, 2019 Andy Aschwanden and Constantine Khroulev
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

#include "DischargeRouting.hh"

#include "pism/util/IceGrid.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh"
#include "FrontalMeltPhysics.hh"

namespace pism {
namespace frontalmelt {
  
DischargeRouting::DischargeRouting(IceGrid::ConstPtr g)
  : CompleteFrontalMeltModel(g, nullptr) {

  m_log->message(2,
             "* Initializing the frontal melt model\n"
             "  UAF-UT\n");
  
  m_frontal_melt_rate = allocate_frontal_melt_rate(g);

  unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");

  m_theta_ocean.reset(new IceModelVec2T(g, "theta_ocean", 1, evaluations_per_year));
}

DischargeRouting::~DischargeRouting() {
  // empty
}

void DischargeRouting::bootstrap_impl(const Geometry &geometry) {
  (void) geometry;

  FrontalMeltPhysics physics(*m_config);

  // FIXME: theta has units of Kelvin, so zero is not a reasonable default value
  m_theta_ocean->set(0.0);
}

void DischargeRouting::init_impl(const Geometry &geometry) {
  (void) geometry;

  FrontalMeltPhysics physics(*m_config);

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.routing");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_theta_ocean = IceModelVec2T::ForcingField(m_grid,
                                                file,
                                                "theta_ocean",
                                                "", // no standard name
                                                buffer_size,
                                                evaluations_per_year,
                                                periodic);
  }

  m_theta_ocean->set_attrs("climate_forcing",
                           "potential temperature of the adjacent ocean",
                           "Celsius", "");

  m_theta_ocean->init(opt.filename, opt.period, opt.reference_time);

}

/*!
 * Initialize potential temperature from IceModelVecs instead of an input
 * file (for testing).
 */
void DischargeRouting::initialize(const IceModelVec2S &theta) {
  m_theta_ocean->copy_from(theta);
}

void DischargeRouting::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {

  m_theta_ocean->update(t, dt);

  FrontalMeltPhysics physics(*m_config);

  const IceModelVec2CellType &cell_type           = inputs.geometry->cell_type;
  const IceModelVec2S        &bed_elevation       = inputs.geometry->bed_elevation;
  const IceModelVec2S        &ice_thickness       = inputs.geometry->ice_thickness;
  const IceModelVec2S        &sea_level_elevation = inputs.geometry->sea_level_elevation;
  const IceModelVec2S        &water_flux          = *inputs.subglacial_water_flux;

  IceModelVec::AccessList list
    {&ice_thickness, &bed_elevation, &cell_type, &sea_level_elevation,
     &water_flux, m_theta_ocean.get(), m_frontal_melt_rate.get()};

  double
    seconds_per_day = 86400,
    grid_spacing    = 0.5 * (m_grid->dx() + m_grid->dy());

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.icy(i, j)) {

      // Assume for now that thermal forcing is equal to theta_ocean. Also, thermal
      // forcing is generally not available at the grounding line.
      //
      double TF = (*m_theta_ocean)(i, j);

      double cross_section_area = ice_thickness(i, j) * grid_spacing;

      // subglacial discharge: convert from m^2/s to m/day
      //
      // The parameterization uses the water "flux" in the very strange units of "m /
      // day". The only interpretation of it I can come up with is that it is the rate of
      // retreat of the front that results in the mass loss equivalent to the mass loss
      // due to the water flux from underneath the ice.
      //
      // [flux] = m^2 / s, so
      // [flux * grid_spacing] = m^3 / s, so
      // [flux * grid_spacing / cross_section_area] = m / s, and
      // [flux * grid_spacing  * (s / day) / cross_section_area] = m / day
      double Q_sg = water_flux(i, j) * grid_spacing;
      double q_sg = Q_sg / cross_section_area * seconds_per_day;

      // get the average ice thickness over ice-covered grounded neighbors
      double water_depth = sea_level_elevation(i, j) - bed_elevation(i, j);

      (*m_frontal_melt_rate)(i, j) = physics.frontal_melt_from_undercutting(water_depth, q_sg, TF);
      // convert from m / day to m / s
      (*m_frontal_melt_rate)(i, j) /= seconds_per_day;
    } else {

      // This parameterization is applicable at grounded termini (see the case above), but
      // *not* at calving fronts of ice shelves.
      (*m_frontal_melt_rate)(i, j) = 0.0;
    }
  } // end of the loop over grid points
}

MaxTimestep DischargeRouting::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("frontalmelt routing");
}


const IceModelVec2S& DischargeRouting::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace frontalmelt
} // end of namespace pism
