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

#ifndef _PSTEMPERATUREINDEX_H_
#define _PSTEMPERATUREINDEX_H_

#include <memory>

#include "pism/util/iceModelVec2T.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "localMassBalance.hh"

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
class TemperatureIndex : public SurfaceModel {
public:
  TemperatureIndex(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> input);
  virtual ~TemperatureIndex();

  // diagnostics (for the last time step)
  const IceModelVec2S& firn_depth() const;
  const IceModelVec2S& snow_depth() const;
  // these represent totals (not rates) over the time step
  const IceModelVec2S& air_temp_sd() const;
  const IceModelVec2S& accumulation() const;
  const IceModelVec2S& melt() const;
  const IceModelVec2S& runoff() const;

protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  virtual const IceModelVec2S& mass_flux_impl() const;
  virtual const IceModelVec2S& temperature_impl() const;

  double compute_next_balance_year_start(double time);
protected:
  //! mass balance scheme to use

  std::unique_ptr<LocalMassBalance> m_mbscheme;

  //! if not NULL then user wanted fausto PDD stuff
  std::unique_ptr<FaustoGrevePDDObject> m_faustogreve;

  //! holds degree-day factors in location-independent case
  LocalMassBalance::DegreeDayFactors m_base_ddf;

  //! K; daily amount of randomness
  double m_base_pddStdDev;

  double m_next_balance_year_start;

  //! cached surface mass balance rate
  IceModelVec2S m_mass_flux;

  IceModelVec2S::Ptr m_temperature;

  //! firn depth
  IceModelVec2S m_firn_depth;

  //! snow depth (reset once a year)
  IceModelVec2S m_snow_depth;

  //! standard deviation of the daily variability of the air temperature
  IceModelVec2T::Ptr m_air_temp_sd;

  //! total accumulation during the last time step
  IceModelVec2S m_accumulation;

  //! total melt during the last time step
  IceModelVec2S m_melt;

  //! total runoff during the last time step
  IceModelVec2S m_runoff;

  bool m_sd_use_param, m_sd_file_set;
  int m_sd_period;
  double m_sd_param_a, m_sd_param_b;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSTEMPERATUREINDEX_H_ */
