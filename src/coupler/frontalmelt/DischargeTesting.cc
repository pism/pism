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

#include "DischargeTesting.hh"

#include "pism/util/IceGrid.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/coupler/util/options.hh"
#include "FrontalMeltPhysics.hh"

namespace pism {
namespace frontalmelt {
  
DischargeTesting::DischargeTesting(IceGrid::ConstPtr g)
  : CompleteFrontalMeltModel(g, nullptr) {

  m_log->message(2,
             "* Initializing the frontal melt model\n"
             "  UAF-UT testing\n");
  
  m_frontal_melt_rate = allocate_frontal_melt_rate(g);

  unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");

  m_theta_ocean.reset(new IceModelVec2T(g, "theta_ocean", 1, evaluations_per_year));
  m_salinity_ocean.reset(new IceModelVec2T(g, "salinity_ocean", 1, evaluations_per_year));
}

DischargeTesting::~DischargeTesting() {
  // empty
}

void DischargeTesting::bootstrap_impl(const Geometry &geometry) {
  (void) geometry;

  FrontalMeltPhysics physics(*m_config);

  m_theta_ocean->set(0.0);
  m_salinity_ocean->set(0.0);
}

void DischargeTesting::init_impl(const Geometry &geometry) {
  (void) geometry;

  FrontalMeltPhysics physics(*m_config);
  
  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.testing");

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
void DischargeTesting::initialize(const IceModelVec2S &theta, const IceModelVec2S &salinity) {
  m_theta_ocean->copy_from(theta);
  m_salinity_ocean->copy_from(salinity);
}

void DischargeTesting::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {
  (void) t;
  (void) dt;

  FrontalMeltPhysics physics(*m_config);

  // We implement Eq. 1 from Rignot et al (2016):
  // q_m = (A * h * Q_sg ^{\alpha} + B) * TF^{\beta}; where
  // A, B, alpha, beta are tuning parameters
  // Rignot (2016) is an update on Xu 2013
  
  double water_density = m_config->get_double("constants.fresh_water.density");

  const IceModelVec2CellType &cell_type = inputs.geometry->cell_type;
  // ice thickness, meters
  const IceModelVec2S &ice_thickness = inputs.geometry->ice_thickness;
  // cell area, meters^2
  const IceModelVec2S &cell_area = inputs.geometry->cell_area;
  // subglacial discharge, mass change over this time step
  const IceModelVec2S &discharge = *inputs.subglacial_discharge;

  IceModelVec::AccessList list
    {&ice_thickness, &cell_type, &cell_area, &discharge, m_theta_ocean.get(),
     m_salinity_ocean.get(), m_frontal_melt_rate.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double h = ice_thickness(i,j);

    // only use grounded cell with non-zero ice thickness are considered
    if (cell_type.grounded(i, j) and (h > 0.0)) {
        
        // Assume for now that thermal forcing is equal to theta_ocean also, thermal forcing
        // is generally not available at the grounding line.
        double TF = (*m_theta_ocean)(i, j);
        
        // subglacial discharge: convert from kg to m/s
        double Qsg = discharge(i, j) / (water_density * cell_area(i, j) * dt);
        
        (*m_frontal_melt_rate)(i,j) = physics.frontal_melt_from_undercutting(h, Qsg, TF);
      } else {

      (*m_frontal_melt_rate)(i,j) = 0.0;

    }
      
  } // end of the loop over grid points
}

MaxTimestep DischargeTesting::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("frontalmelt testing");
}


const IceModelVec2S& DischargeTesting::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace frontalmelt
} // end of namespace pism
