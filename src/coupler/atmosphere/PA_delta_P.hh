// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#ifndef _PADPFORCING_H_
#define _PADPFORCING_H_

#include "PScalarForcing.hh"
#include "PAModifier.hh"

namespace pism {

class PA_delta_P : public PScalarForcing<AtmosphereModel,PAModifier>
{
public:
  PA_delta_P(IceGrid &g, const Config &conf, AtmosphereModel* in);
  virtual ~PA_delta_P();

  virtual void init(Vars &vars);
  virtual void init_timeseries(const std::vector<double> &ts);

  virtual void mean_precipitation(IceModelVec2S &result);

  virtual void precip_time_series(int i, int j, std::vector<double> &values);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);

  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);

protected:
  NCSpatialVariable air_temp, precipitation;
  std::vector<double> m_offset_values;
private:
  PetscErrorCode allocate_PA_delta_P();
};

} // end of namespace pism

#endif /* _PADPFORCING_H_ */
