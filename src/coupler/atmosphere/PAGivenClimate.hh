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

#ifndef _PAGIVEN_H_
#define _PAGIVEN_H_

#include "PAModifier.hh"
#include "PGivenClimate.hh"

namespace pism {

class PAGivenClimate : public PGivenClimate<PAModifier,AtmosphereModel>
{
public:
  PAGivenClimate(const IceGrid &g);
  virtual ~PAGivenClimate();

  virtual void init();
  virtual void update(double my_t, double my_dt);

  virtual void mean_precipitation(IceModelVec2S &result);
  virtual void mean_annual_temp(IceModelVec2S &result); 
  virtual void temp_snapshot(IceModelVec2S &result);

  virtual void begin_pointwise_access();
  virtual void end_pointwise_access();

  virtual void init_timeseries(const std::vector<double> &ts);
  virtual void temp_time_series(int i, int j, std::vector<double> &values);
  virtual void precip_time_series(int i, int j, std::vector<double> &values);
protected:
  IceModelVec2T *precipitation, *air_temp;
};

} // end of namespace pism

#endif /* _PAGIVEN_H_ */
