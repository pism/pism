// Copyright (C) 2008-2018 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PASeariseGreenland_hh
#define __PASeariseGreenland_hh

#include "YearlyCycle.hh"
#include "pism/util/Timeseries.hh"

namespace pism {
namespace atmosphere {

//! \brief A modification of YearlyCycle tailored for the
//! SeaRISE-Greenland assessment. Uses the Fausto [\ref Faustoetal2009]
//! present-day temperature parameterization and stored precipitation data.
class SeaRISEGreenland : public YearlyCycle {
public:
  SeaRISEGreenland(IceGrid::ConstPtr g);
  virtual ~SeaRISEGreenland();

  virtual void init_impl(const Geometry &geometry);
  virtual void precip_time_series_impl(int i, int j, std::vector<double> &values) const;
protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);
};


} // end of namespace atmosphere
} // end of namespace pism

#endif  // __PASeariseGreenland_hh
