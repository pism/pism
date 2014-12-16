// Copyright (C) 2011, 2013, 2014 Constantine Khroulev
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

#ifndef _PODIRECTFORCING_H_
#define _PODIRECTFORCING_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"

namespace pism {

class POGiven : public PGivenClimate<POModifier,OceanModel>
{
public:
  POGiven(IceGrid &g);
  virtual ~POGiven();

  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);

  virtual void sea_level_elevation(double &result);

  virtual void shelf_base_temperature(IceModelVec2S &result);
  virtual void shelf_base_mass_flux(IceModelVec2S &result);
  virtual void melange_back_pressure_fraction(IceModelVec2S &result);
protected:
  IceModelVec2T *shelfbtemp, *shelfbmassflux;
};


} // end of namespace pism

#endif /* _PODIRECTFORCING_H_ */
