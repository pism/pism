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

#ifndef __PISMAtmosphere_hh
#define __PISMAtmosphere_hh

#include "PISMComponent.hh"

namespace pism {
class IceModelVec2S;

///// AtmosphereModel: models which provide precipitation and temperature
/////                      to the SurfaceModel below

//! A purely virtual class defining the interface of a PISM Atmosphere Model.
class AtmosphereModel : public Component_TS {
public:
  AtmosphereModel(IceGrid &g, const Config &conf)
    : Component_TS(g, conf) {};

  //! \brief Sets result to the mean precipitation, in m/s ice equivalent.
  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result) = 0;

  //! \brief Sets result to the mean annual near-surface air temperature, in degrees Kelvin.
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result) = 0;

  virtual PetscErrorCode begin_pointwise_access() = 0;
  virtual PetscErrorCode end_pointwise_access() = 0;
  virtual PetscErrorCode init_timeseries(const std::vector<double> &ts) = 0;
  //! \brief Sets a pre-allocated N-element array "result" to the time-series of
  //! ice-equivalent precipitation (m/s) at the point i,j on the grid.
  //!
  //! See temp_time_series() for more.
  virtual PetscErrorCode precip_time_series(int i, int j, std::vector<double> &result) = 0;

  //! \brief Sets a pre-allocated N-element array "result" to the time-series
  //! of near-surface air temperature (degrees Kelvin) at the point i,j on the
  //! grid. Times (in years) are specified in ts. NB! Has to be surrounded by
  //! begin_pointwise_access() and end_pointwise_access()
  virtual PetscErrorCode temp_time_series(int i, int j, std::vector<double> &result) = 0;
  //! \brief Sets result to a snapshot of temperature for the current time.
  //! (For diagnostic purposes.)
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result) = 0;
protected:

  std::vector<double> m_ts_times;
};

} // end of namespace pism

#endif  // __PISMAtmosphere_hh
