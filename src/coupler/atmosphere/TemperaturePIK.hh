// Copyright (C)  2009-2018 Ricarda Winkelmann, Torsten Albrecht

//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef __PATemperaturePIK_hh
#define __PATemperaturePIK_hh

#include "YearlyCycle.hh"
#include "pism/util/Timeseries.hh"

namespace pism {
namespace atmosphere {


class TemperaturePIK : public YearlyCycle 
{
public:
  TemperaturePIK(IceGrid::ConstPtr g);
  virtual ~TemperaturePIK();

  virtual void init_impl(const Geometry &geometry);
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &values) const;

protected:

  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);
 
  bool temp_huybrechts_dewolde99_set, temp_era_interim_set, temp_era_interim_sin_set, temp_era_interim_lon_set;

};

} // end of namespace atmosphere
} // end of namespace pism

#endif	// __PATemperaturePIK_hh
