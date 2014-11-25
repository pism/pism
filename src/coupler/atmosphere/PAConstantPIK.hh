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

#ifndef _PACONSTANTPIK_H_
#define _PACONSTANTPIK_H_

#include "iceModelVec.hh"
#include "PISMAtmosphere.hh"

namespace pism {

class PAConstantPIK : public AtmosphereModel {
public:
  PAConstantPIK(IceGrid &g, const Config &conf);
  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);
  virtual void mean_precipitation(IceModelVec2S &result);
  virtual void mean_annual_temp(IceModelVec2S &result);
  virtual void begin_pointwise_access();
  virtual void end_pointwise_access();
  virtual void temp_time_series(int i, int j, std::vector<double> &values);
  virtual void precip_time_series(int i, int j, std::vector<double> &values);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
  virtual void temp_snapshot(IceModelVec2S &result);
  virtual void init_timeseries(const std::vector<double> &ts);
protected:
  IceModelVec2S *usurf, *lat;
  std::string input_file;
  IceModelVec2S precipitation, air_temp;
  NCSpatialVariable air_temp_snapshot;
};

} // end of namespace pism

#endif /* _PACONSTANTPIK_H_ */
