// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include "base/util/iceModelVec2T.hh"
#include "coupler/PISMSurface.hh"
#include "localMassBalance.hh"

namespace pism {
namespace surface {

//! \brief A class implementing a temperature-index (positive degree-day) scheme
//! to compute melt and runoff, and thus surface mass balance, from
//! precipitation and air temperature.
/*! 
  Temperature-index schemes are far from perfect as a way of modeling surface mass
  balance on ice sheets which experience surface melt, but they are known to have
  reasonable data requirements and to do a good job when tuned appropriately
  [\ref Hock05].

  This base class already accesses a fair amount of functionality.  It holds a
  pointer to an instance of the LocalMassBalance class.  This class has method
  LocalMassBalance::getMassFluxFromTemperatureTimeSeries() which uses the
  precipitation during the ice sheet model time step, plus a variable temperature
  over that time step, to compute melt, refreeze, and surface balance.
*/
class TemperatureIndex : public SurfaceModel {
public:
  TemperatureIndex(IceGrid::ConstPtr g);
  virtual ~TemperatureIndex();

  const IceModelVec2S& surface_accumulation() const;
  const IceModelVec2S& surface_melt() const;
  const IceModelVec2S& surface_runoff() const;
  const IceModelVec2S& snow_depth() const;
  const IceModelVec2S& air_temp_sd() const;

  const IceModelVec2S& cumulative_surface_accumulation() const;
  const IceModelVec2S& cumulative_surface_melt() const;
  const IceModelVec2S& cumulative_surface_runoff() const;

protected:
  virtual void init_impl();
  virtual void update_impl(double my_t, double my_dt);

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result) const;
  virtual void ice_surface_temperature_impl(IceModelVec2S &result) const;
  virtual MaxTimestep max_timestep_impl(double t) const;

  double compute_next_balance_year_start(double time);
protected:
  LocalMassBalance *m_mbscheme;         //!< mass balance scheme to use

  FaustoGrevePDDObject *m_faustogreve;  //!< if not NULL then user wanted fausto PDD stuff

  //! holds degree-day factors in location-independent case
  LocalMassBalance::DegreeDayFactors m_base_ddf;

  //! K; daily amount of randomness
  double m_base_pddStdDev;
  //! K; temps are positive above this
  double m_base_pddThresholdTemp;
  double m_next_balance_year_start;

  //! cached surface mass balance rate
  IceModelVec2S m_climatic_mass_balance;

  //! diagnostic output accumulation rate (snow - rain)
  IceModelVec2S m_accumulation_rate;

  //! diagnostic output melt rate (rate at which snow and ice is melted; but some snow melt
  //! refreezes)
  IceModelVec2S m_melt_rate;

  //! diagnostic output meltwater runoff rate
  IceModelVec2S m_runoff_rate;

  //! snow depth (reset once a year)
  IceModelVec2S m_snow_depth;
  IceModelVec2T m_air_temp_sd;

  IceModelVec2S m_cumulative_accumulation;
  IceModelVec2S m_cumulative_melt;
  IceModelVec2S m_cumulative_runoff;

  bool m_randomized, m_randomized_repeatable, m_use_fausto_params;
  bool m_sd_use_param, m_sd_file_set;
  int m_sd_period, m_sd_period_years;
  double m_sd_ref_time, m_sd_param_a, m_sd_param_b;
};

/*! @brief Surface accumulation rate. */
class PDD_saccum : public Diag<TemperatureIndex>
{
public:
  PDD_saccum(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Surface melt rate. */
class PDD_smelt : public Diag<TemperatureIndex>
{
public:
  PDD_smelt(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Surface runoff rate. */
class PDD_srunoff : public Diag<TemperatureIndex>
{
public:
  PDD_srunoff(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Snow cover depth. */
class PDD_snow_depth : public Diag<TemperatureIndex>
{
public:
  PDD_snow_depth(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl();
};

/*! @brief Standard deviation of near-surface air temperature. */
class PDD_air_temp_sd : public Diag<TemperatureIndex>
{
public:
  PDD_air_temp_sd(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl();
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSTEMPERATUREINDEX_H_ */
