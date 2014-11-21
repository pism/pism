// Copyright (C) 2008-2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "PAYearlyCycle.hh"
#include "Timeseries.hh"

namespace pism {

//! \brief A modification of PAYearlyCycle tailored for the
//! SeaRISE-Greenland assessment. Uses the Fausto [\ref Faustoetal2009]
//! present-day temperature parameterization and stored precipitation data.
class PA_SeaRISE_Greenland : public PAYearlyCycle {
public:
  PA_SeaRISE_Greenland(IceGrid &g, const Config &conf);
  virtual ~PA_SeaRISE_Greenland();

  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);
  virtual void precip_time_series(int i, int j, std::vector<double> &values);
protected:
  IceModelVec2S *m_lat, *m_lon, *m_surfelev;
};


} // end of namespace pism

#endif  // __PASeariseGreenland_hh
