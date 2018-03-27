// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PAGENERICYEARLYCYCLE_H_
#define _PAGENERICYEARLYCYCLE_H_

#include <memory>               // unique_ptr

#include "YearlyCycle.hh"

namespace pism {
class Timeseries;

namespace atmosphere {

class CosineYearlyCycle : public YearlyCycle {
public:
  CosineYearlyCycle(IceGrid::ConstPtr g);
  virtual ~CosineYearlyCycle();

  virtual void init_impl(const Geometry &geometry);
  virtual void init_timeseries_impl(const std::vector<double> &ts) const;
protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  std::unique_ptr<Timeseries> m_A;                 // amplitude scaling
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAGENERICYEARLYCYCLE_H_ */
