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

#ifndef _PAANOMALY_H_
#define _PAANOMALY_H_

#include "PGivenClimate.hh"
#include "PAModifier.hh"

namespace pism {

//! \brief Reads and uses air_temp and precipitation anomalies from a file.
class PAAnomaly : public PGivenClimate<PAModifier,AtmosphereModel>
{
public:
  PAAnomaly(IceGrid &g, const Config &conf, AtmosphereModel* in);
  virtual ~PAAnomaly();

  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);

  virtual void mean_precipitation(IceModelVec2S &result);
  virtual void mean_annual_temp(IceModelVec2S &result); 
  virtual void temp_snapshot(IceModelVec2S &result);

  virtual void init_timeseries(const std::vector<double> &ts);
  virtual void begin_pointwise_access();
  virtual void end_pointwise_access();
  virtual void temp_time_series(int i, int j, std::vector<double> &values);
  virtual void precip_time_series(int i, int j, std::vector<double> &values);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);

  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);

protected:
  std::vector<double> ts_mod, ts_values;
  NCSpatialVariable air_temp, precipitation;
  IceModelVec2T *air_temp_anomaly, *precipitation_anomaly;
  std::vector<double> m_mass_flux_anomaly, m_temp_anomaly;
private:
  PetscErrorCode allocate_PAAnomaly();
};

} // end of namespace pism

#endif /* _PAANOMALY_H_ */
