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

#ifndef _PALAPSERATES_H_
#define _PALAPSERATES_H_

#include "coupler/util/PLapseRates.hh"
#include "PAModifier.hh"

namespace pism {
namespace atmosphere {

class LapseRates : public PLapseRates<AtmosphereModel,PAModifier>
{
public:
  LapseRates(IceGrid::ConstPtr g, AtmosphereModel* in);
  virtual ~LapseRates();

protected:
  virtual void init_impl();

  virtual void mean_precipitation_impl(IceModelVec2S &result) const;
  virtual void mean_annual_temp_impl(IceModelVec2S &result) const;

  virtual void begin_pointwise_access_impl() const;
  virtual void end_pointwise_access_impl() const;

  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &result) const;
  virtual void temp_time_series_impl(int i, int j, std::vector<double> &result) const;

protected:
  double m_precip_lapse_rate;
  const IceModelVec2S *m_surface;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PALAPSERATES_H_ */

