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

#ifndef _PAMODIFIER_H_
#define _PAMODIFIER_H_

#include "PISMAtmosphere.hh"

namespace pism {

class PAModifier : public Modifier<AtmosphereModel>
{
public:
  PAModifier(const IceGrid &g, AtmosphereModel* in)
    : Modifier<AtmosphereModel>(g, in) {}
  virtual ~PAModifier() {}

  virtual void mean_precipitation(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->mean_precipitation(result);
    }
  }

  virtual void mean_annual_temp(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->mean_annual_temp(result);
    }
  }

  virtual void begin_pointwise_access()
  {
    if (input_model != NULL) {
      input_model->begin_pointwise_access();
    }
  }

  virtual void end_pointwise_access()
  {
    if (input_model != NULL) {
      input_model->end_pointwise_access();
    }
  }

  virtual void temp_time_series(int i, int j, std::vector<double> &result)
  {
    if (input_model != NULL) {
      input_model->temp_time_series(i, j, result);
    }
  }

  virtual void precip_time_series(int i, int j, std::vector<double> &result)
  {
    if (input_model != NULL) {
      input_model->precip_time_series(i, j, result);
    }
  }

  virtual void temp_snapshot(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->temp_snapshot(result);
    }
  }

  virtual void init_timeseries(const std::vector<double> &ts)
  {
    if (input_model != NULL) {
      input_model->init_timeseries(ts);
    }

    m_ts_times = ts;
  }
};

} // end of namespace pism

#endif /* _PAMODIFIER_H_ */
