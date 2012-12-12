// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _PSSIMPLE_H_
#define _PSSIMPLE_H_

#include "PISMSurface.hh"
#include "PISMAtmosphere.hh"
#include "NCSpatialVariable.hh"

//! \brief A class implementing a primitive surface model.
/*! 
This is an "invisible" surface processes model which "passes through"
information from the atmosphere above directly to the ice below the surface
layers.  It implements two modeling choices:
  \li accumulation which is obtained from an atmosphere model is interpreted
      as surface mass flux;
  \li mean-annual near-surface air temperature is interpreted as instantaneous
      temperature of the ice at the ice surface.

The second choice means that the upper boundary condition of the conservation of
energy scheme for the ice fluid is exactly the 2m air temperature.
*/
class PSSimple : public PISMSurfaceModel {
public:
  PSSimple(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf) {};
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  {
    t = my_t; dt = my_dt;
    if (atmosphere) {
      PetscErrorCode ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);
    }
    return 0;
  }
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);
protected:
  NCSpatialVariable climatic_mass_balance, ice_surface_temp;
};

#endif /* _PSSIMPLE_H_ */
