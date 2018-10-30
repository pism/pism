// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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
  m_salinity_ocean.reset(new IceModelVec2T(g, "salinity_ocean", 1, evaluations_per_year));
}

DischargeRouting::~DischargeRouting() {
  // empty
}

void DischargeRouting::bootstrap_impl(const Geometry &geometry) {
  (void) geometry;

  FrontalMeltPhysics physics(*m_config);

  m_theta_ocean->set(0.0);
  m_salinity_ocean->set(0.0);
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

    m_salinity_ocean = IceModelVec2T::ForcingField(m_grid,
                                                   file,
                                                   "salinity_ocean",
                                                   "", // no standard name
                                                   buffer_size,
                                                   evaluations_per_year,
                                                   periodic);
  }

  m_theta_ocean->set_attrs("climate_forcing",
                           "potential temperature of the adjacent ocean",
                           "Kelvin", "");

  m_salinity_ocean->set_attrs("climate_forcing",
                              "salinity of the adjacent ocean",
                              "g/kg", "");
}

/*!
 * Initialize potential temperature and salinity from IceModelVecs instead of an input
 * file (for testing).
 */
void DischargeRouting::initialize(const IceModelVec2S &theta, const IceModelVec2S &salinity) {
  m_theta_ocean->copy_from(theta);
  m_salinity_ocean->copy_from(salinity);
}

void DischargeRouting::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {
  (void) t;
  (void) dt;

  FrontalMeltPhysics physics(*m_config);
    
  double water_density = m_config->get_double("constants.fresh_water.density");

  const IceModelVec2CellType &cell_type = inputs.geometry->cell_type;
  // ice thickness, meters
  const IceModelVec2S &ice_thickness = inputs.geometry->ice_thickness;
  // cell area, meters^2
  const IceModelVec2S &cell_area = inputs.geometry->cell_area;
  // subglacial discharge, mass change over this time step
  const IceModelVec2S &subglacial_discharge = *inputs.subglacial_discharge_at_grounding_line;

  IceModelVec::AccessList list
    {&ice_thickness, &cell_type, &cell_area, &subglacial_discharge, m_theta_ocean.get(),
     m_salinity_ocean.get(), m_frontal_melt_rate.get()};

  // index offsets for iterating over neighbors
  const int i_offsets[4] = {1, 0, -1, 0};
  const int j_offsets[4] = {0, 1, 0, -1};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ocean(i, j) and cell_type.next_to_grounded_ice(i, j)) {

      // Assume for now that thermal forcing is equal to theta_ocean also, thermal forcing
      // is generally not available at the grounding line.
      double TF = (*m_theta_ocean)(i, j);

      // subglacial discharge: convert from kg to m/s
      double Qsg = subglacial_discharge(i, j) / (water_density * cell_area(i, j) * dt);

      // get the average ice thickness over ice-covered grounded neighbors
      double H = 0.0;
      {
        int n_grounded_neighbors = 0;
        for (int k = 0; k < 4; ++k) {
          int i_n = i + i_offsets[k];
          int j_n = j + j_offsets[k];

          if (cell_type.grounded_ice(i_n, j_n)) {
            H += ice_thickness(i_n, j_n);
            n_grounded_neighbors += 1;
          }
        }

        if (n_grounded_neighbors > 0) {
          H /= n_grounded_neighbors;
        }
      }

      (*m_frontal_melt_rate)(i,j) = physics.frontal_melt_from_undercutting(H, Qsg, TF);
    } else { // end of "if this is an ocean cell next to grounded ice"

      // This parameterization is applicable at grounded termini (see the case above), but
      // *not* at calving fronts of ice shelves.
      (*m_frontal_melt_rate)(i,j) = 0.0;
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
