// Copyright (C) 2008-2011, 2013, 2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#ifndef __PISMOceanModel_hh
#define __PISMOceanModel_hh

#include "PISMComponent.hh"

namespace pism {

class IceModelVec2S;

//! A very rudimentary PISM ocean model.
class PISMOceanModel : public PISMComponent_TS {
public:
  PISMOceanModel(IceGrid &g, const PISMConfig &conf)
    : PISMComponent_TS(g, conf), sea_level(0) {}
  virtual ~PISMOceanModel() {};

  virtual PetscErrorCode sea_level_elevation(double &result) = 0;
  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result) = 0;
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result) = 0;
  virtual PetscErrorCode melange_back_pressure_fraction(IceModelVec2S &result);
protected:
  double sea_level;
};

} // end of namespace pism

#endif  // __PISMOceanModel_hh
