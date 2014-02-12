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

//! \brief Reads and uses air_temp and precipitation anomalies from a file.
class PAAnomaly : public PGivenClimate<PAModifier,PISMAtmosphereModel>
{
public:
  PAAnomaly(IceGrid &g, const PISMConfig &conf, PISMAtmosphereModel* in);
  virtual ~PAAnomaly();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result); 
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual PetscErrorCode init_timeseries(double *ts, unsigned int N);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, double *values);
  virtual PetscErrorCode precip_time_series(int i, int j, double *values);

  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);

protected:
  std::vector<double> ts_mod, ts_values;
  NCSpatialVariable air_temp, precipitation;
  IceModelVec2T *air_temp_anomaly, *precipitation_anomaly;
  std::vector<double> m_mass_flux_anomaly, m_temp_anomaly;
private:
  PetscErrorCode allocate_PAAnomaly();
};

#endif /* _PAANOMALY_H_ */
