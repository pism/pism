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

#ifndef _PAGENERICYEARLYCYCLE_H_
#define _PAGENERICYEARLYCYCLE_H_

#include "PAYearlyCycle.hh"

namespace pism {

class Timeseries;

class PACosineYearlyCycle : public PAYearlyCycle {
public:
  PACosineYearlyCycle(IceGrid &g, const Config &conf);
  virtual ~PACosineYearlyCycle();

  virtual void init(Vars &vars);
  virtual void init_timeseries(const std::vector<double> &ts);
  virtual void update(double my_t, double my_dt);
  virtual void temp_snapshot(IceModelVec2S &result);
protected:
  Timeseries *A;                 // amplitude scaling
};

} // end of namespace pism

#endif /* _PAGENERICYEARLYCYCLE_H_ */
