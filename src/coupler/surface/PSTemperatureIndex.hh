// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "iceModelVec2T.hh"
#include "PISMSurface.hh"
#include "localMassBalance.hh"
#include "NCVariable.hh"

namespace pism {

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
class PSTemperatureIndex : public SurfaceModel {
public:
  PSTemperatureIndex(IceGrid &g);
  virtual ~PSTemperatureIndex();
  virtual void update(double my_t, double my_dt);
  virtual void init();
  virtual void max_timestep(double my_t, double &my_dt, bool &restrict);
  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);  
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  LocalMassBalance *mbscheme;         //!< mass balance scheme to use

  FaustoGrevePDDObject *faustogreve;  //!< if not NULL then user wanted fausto PDD stuff

  DegreeDayFactors base_ddf;          //!< holds degree-day factors in location-independent case
  double  base_pddStdDev,        //!< K; daily amount of randomness
    base_pddThresholdTemp, //!< K; temps are positive above this
    m_next_balance_year_start;
  IceModelVec2S
  climatic_mass_balance, //!< cached surface mass balance rate
    accumulation_rate,     //!< diagnostic output accumulation rate (snow - rain)
    melt_rate,             //!< diagnostic output melt rate (rate at which snow
  //!< and ice is melted, but some snow melt refreezes)
    runoff_rate,           //!< diagnostic output meltwater runoff rate
    snow_depth;            //!< snow depth (reset once a year)
  IceModelVec2T air_temp_sd;

  IceModelVec2S *lat, *lon, *usurf;
  //!< PSTemperatureIndex must hold these pointers in order to use
  //! object which needs 3D location to determine degree day factors.
  IceModelVec2Int *mask;

  NCSpatialVariable ice_surface_temp;

  bool randomized, randomized_repeatable, fausto_params;
  bool sd_file_set, sd_period_set, sd_ref_year_set, sd_use_param;
  int sd_period, sd_period_years, sd_ref_year;
  double sd_ref_time, sd_param_a, sd_param_b;
  std::string filename;
  double compute_next_balance_year_start(double time);
};

} // end of namespace pism

#endif /* _PSTEMPERATUREINDEX_H_ */
