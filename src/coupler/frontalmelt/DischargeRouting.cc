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
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"

#include "pism/coupler/util/options.hh"

namespace pism {
namespace frontalmelt {
  
DischargeRouting::DischargeRouting(IceGrid::ConstPtr g)
  : FrontalMeltModel(g, nullptr) {

  m_frontal_melt_rate   = allocate_frontal_melt_rate(g);

}

DischargeRouting::~DischargeRouting() {
  // empty
}

void DischargeRouting::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
             "* Initializing the frontal melt model UAF\n"
             "  from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.routing");

  // This model requires hydrology "routing" or "distributed"

  std::string method = m_config->get_string("hydrology.model");

  if (method != "routing" || method != "distributed") {
    throw RuntimeError(PISM_ERROR_LOCATION, "option -frontal_melt routing"
                       " requires -hydrology {routing,distributed}");
  }
  
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


  // read time-independent data right away:
  if (m_theta_ocean->n_records() == 1 && m_salinity_ocean->n_records() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void DischargeRouting::bootstrap_impl(const Geometry &geometry) {
  (void) geometry;
}
  
void DischargeRouting::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;
  (void) geometry;

  // We implement Eq. 1 from Rignot et al (2016):
  // q_m = (A * h * Q_sg ^{\alpha} + B) * TF^{\beta}; where
  // A, B, alpha, beta are tuning parameters
  // Rignot (2016) is an update on Xu 2013
  
  double A = m_config->get_double("frontal_melt.parameter_a");
  double B = m_config->get_double("frontal_melt.parameter_b");
  double alpha = m_config->get_double("frontal_melt.power_alpha");
  double beta = m_config->get_double("frontal_melt.power_beta");
  
  // get ice thickness
  const IceModelVec2S* ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  // get subglacial discharge
  const IceModelVec2S *subglacial_discharge  = m_grid->variables().get_2d_scalar("tendency_of_subglacial_water_mass_at_grounding_line");

  IceModelVec::AccessList list{ ice_thickness, subglacial_discharge };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double h = (*ice_thickness)(i,j);
    // Assume for now that thermal forcing is equal to theta_coean
    // also, thermal forcing is generally not available at the grounding line.
    double TF = (*m_theta_ocean)(i,j);
    double Qsg = (*subglacial_discharge)(i, j);
    (*m_frontal_melt_rate)(i,j) = (A * h * pow(Qsg, alpha) + B) * pow(TF, beta);
  }

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
