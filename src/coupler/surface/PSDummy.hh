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

#ifndef _PSDUMMY_H_
#define _PSDUMMY_H_

#include "PISMSurface.hh"
#include "PISMAtmosphere.hh"

//! A do-nothing (dummy) surface model. <b> Please avoid using it for real modeling! </b>
/*!
  This dummy class is used, for example, when an internal (with respect to
  IceModel) formula generates the surface mass balance. A specific case is the
  manufactured solutions used in verification.
*/
class PSDummy : public PISMSurfaceModel {
public:
  PSDummy(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {};

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode init(PISMVars &) { return 0; };
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  { t = my_t; dt = my_dt; return 0; } // do nothing
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S&)
  { return 0; }

  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &)
  { return 0; }
  virtual void add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &/*result*/) {}
  virtual PetscErrorCode define_variables(set<string> /*vars*/, const PIO &/*nc*/, PISM_IO_Type /*nctype*/)
  { return 0; }
  virtual PetscErrorCode write_variables(set<string>, const PIO &)
  { return 0; }

  // Does not have an atmosphere model.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}
};

#endif /* _PSDUMMY_H_ */
