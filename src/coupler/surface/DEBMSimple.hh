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
#include "pism/util/ScalarForcing.hh"
#include "pism/util/iceModelVec2T.hh"

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
  const IceModelVec2S &firn_depth() const;
  const IceModelVec2S &snow_depth() const;
  // these represent totals (not rates) over the time step
  const IceModelVec2S &air_temp_sd() const;
  const IceModelVec2S &accumulation_impl() const;
  const IceModelVec2S &melt_impl() const;
  const IceModelVec2S &runoff_impl() const;
  //???
  const IceModelVec2S &surface_insolation_melt() const;
  const IceModelVec2S &surface_temperature_melt() const;
  const IceModelVec2S &surface_offset_melt() const;
  const IceModelVec2S &albedo() const;
  const IceModelVec2S &transmissivity() const;
  const IceModelVec2S &TOAinsol() const;
  const IceModelVec2S &qinsol() const;

private:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const File &output) const;
  virtual void write_model_state_impl(const File &output) const;

  virtual DiagnosticList diagnostics_impl() const;

  virtual const IceModelVec2S &mass_flux_impl() const;
  virtual const IceModelVec2S &temperature_impl() const;

  double compute_next_balance_year_start(double time);

  bool albedo_anomaly_true(double time);

  double earch_sun_distance(double time);
  double earth_declination(double time);

  double earch_sun_distance_paleo(double time);
  double lambda_paleo(double time);
  double earth_declination_paleo(double time);

  DEBMSimplePointwise m_model;

  double m_next_balance_year_start;

  //! cached surface mass balance rate
  IceModelVec2S m_mass_flux;

  IceModelVec2S::Ptr m_temperature;

  //! firn depth
  IceModelVec2S m_firn_depth;

  //! snow depth (reset once a year)
  IceModelVec2S m_snow_depth;

  //! standard deviation of the daily variability of the air temperature
  std::shared_ptr<IceModelVec2T> m_air_temp_sd;

  //! total accumulation during the last time step
  IceModelVec2S::Ptr m_accumulation;

  //! total melt during the last time step
  IceModelVec2S::Ptr m_melt;

  //! total runoff during the last time step
  IceModelVec2S::Ptr m_runoff;

  //! total temperature melt during the last time step
  IceModelVec2S m_tempmelt;

  //! total insolation melt during the last time step
  IceModelVec2S m_insolmelt;

  //! total cmelt during the last timestep
  IceModelVec2S m_cmelt;

  //! albedo field
  IceModelVec2S m_albedo;

  //! if albedo is given as input field
  std::shared_ptr<IceModelVec2T> m_input_albedo;

  //! transmissivity field
  IceModelVec2S m_transmissivity;

  //! TOA insol field
  IceModelVec2S m_TOAinsol;

  //! q insol field
  IceModelVec2S m_qinsol;

  bool m_sd_use_param, m_use_air_temp_sd_file;
  double m_sd_param_a, m_sd_param_b;

  double m_constant_eccentricity;
  double m_constant_perihelion_longitude;
  double m_constant_obliquity;

  std::unique_ptr<ScalarForcing> m_eccentricity;

  std::unique_ptr<ScalarForcing> m_obliquity;

  std::unique_ptr<ScalarForcing> m_perihelion_longitude;

  bool m_use_paleo_file;
};

} // end of namespace surface
} // end of namespace pism

#endif /* PISM_DEBM_SIMPLE_H */
