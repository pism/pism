// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2022 PISM Authors
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

#ifndef PISM_DEBM_SIMPLE_H
#define PISM_DEBM_SIMPLE_H

#include <memory>

#include "DEBMSimplePointwise.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace surface {

//! @brief A class implementing a temperature-index (positive degree-day) scheme
//! to compute melt and runoff, and thus surface mass balance, from
//! precipitation and air temperature.
/*!
  Temperature-index schemes are far from perfect as a way of modeling surface mass
  balance on ice sheets which experience surface melt, but they are known to have
  reasonable data requirements and to do a good job when tuned appropriately
  [@ref Hock05].
*/
class DEBMSimple : public SurfaceModel {
public:
  DEBMSimple(IceGrid::ConstPtr g,
             std::shared_ptr<atmosphere::AtmosphereModel> input);
  virtual ~DEBMSimple() = default;

  // diagnostics (for the last time step)
  const array::Scalar &firn_depth() const;
  const array::Scalar &snow_depth() const;
  // these represent totals (not rates) over the time step
  const array::Scalar &air_temp_sd() const;
  const array::Scalar &accumulation_impl() const;
  const array::Scalar &melt_impl() const;
  const array::Scalar &runoff_impl() const;

  // Contributions to melt from insolation, temperature, and the background melt:
  const array::Scalar &insolation_driven_melt() const;
  const array::Scalar &temperature_driven_melt() const;
  const array::Scalar &background_melt() const;

  // diagnostics
  const array::Scalar &surface_albedo() const;
  const array::Scalar &atmosphere_transmissivity() const;
  const array::Scalar &insolation() const;

private:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;

  virtual const array::Scalar &mass_flux_impl() const;
  virtual const array::Scalar &temperature_impl() const;

  double compute_next_balance_year_start(double time);

  double snow_accumulation(double T, double P) const;

  unsigned int timeseries_length(double dt) const;

  DEBMSimplePointwise m_model;

  double m_next_balance_year_start;

  //! cached surface mass balance rate
  array::Scalar m_mass_flux;

  array::Scalar::Ptr m_temperature;

  //! firn depth
  array::Scalar m_firn_depth;

  //! snow depth (reset once a year)
  array::Scalar m_snow_depth;

  //! standard deviation of the daily variability of the air temperature
  std::shared_ptr<array::Forcing> m_air_temp_sd;

  //! total accumulation during the last time step
  array::Scalar::Ptr m_accumulation;

  //! total melt during the last time step
  array::Scalar::Ptr m_melt;

  //! total runoff during the last time step
  array::Scalar::Ptr m_runoff;

  //! total temperature melt during the last time step
  array::Scalar m_temperature_driven_melt;

  //! total insolation melt during the last time step
  array::Scalar m_insolation_driven_melt;

  //! total background_melt during the last timestep
  array::Scalar m_background_melt;

  //! albedo field
  array::Scalar m_surface_albedo;

  //! if albedo is given as input field
  std::shared_ptr<array::Forcing> m_input_albedo;

  //! transmissivity field
  array::Scalar m_transmissivity;

  //! insolation at the top of the atmosphere
  array::Scalar m_insolation;

  //! year length used to compute the time series length required to get m_n_per_year
  //! evaluations
  double m_year_length;

  //! number of small time steps per year
  unsigned int m_n_per_year;

  bool m_sd_use_param;
  double m_sd_param_a;
  double m_sd_param_b;

  //! interpret all the precipitation as snow (no rain)
  bool m_precip_as_snow;
  //! the temperature below which all precipitation is snow
  double m_Tmin;
  //! the temperature above which all precipitation is rain
  double m_Tmax;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_DEBM_SIMPLE_H */
