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
  PAGivenClimate(IceGrid &g, const Config &conf);
  virtual ~PAGivenClimate();

  virtual PetscErrorCode init(Vars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result); 
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();

  virtual PetscErrorCode init_timeseries(const std::vector<double> &ts);
  virtual PetscErrorCode temp_time_series(int i, int j, std::vector<double> &values);
  virtual PetscErrorCode precip_time_series(int i, int j, std::vector<double> &values);
protected:
  IceModelVec2T *precipitation, *air_temp;
private:
  virtual PetscErrorCode allocate_PAGivenClimate();
};

} // end of namespace pism

#endif /* _PAGIVEN_H_ */
