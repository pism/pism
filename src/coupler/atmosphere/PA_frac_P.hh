// Copyright (C) 2012, 2013, 2014, 2015 PISM Authors
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

#ifndef _PAFPFORCING_H_
#define _PAFPFORCING_H_

#include "PScalarForcing.hh"
#include "PAModifier.hh"

namespace pism {
namespace atmosphere {

class Frac_P : public PScalarForcing<AtmosphereModel,PAModifier>
{
public:
  Frac_P(const IceGrid &g, AtmosphereModel* in);
  virtual ~Frac_P();

  virtual void init();
  virtual void init_timeseries(const std::vector<double> &ts);

  virtual void mean_precipitation(IceModelVec2S &result);

  virtual void precip_time_series(int i, int j, std::vector<double> &values);

protected:
  virtual MaxTimestep max_timestep_impl(double t);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
protected:
  SpatialVariableMetadata air_temp, precipitation;
  std::vector<double> m_offset_values;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAFPFORCING_H_ */
