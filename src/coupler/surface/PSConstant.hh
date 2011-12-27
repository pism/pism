// Copyright (C) 2011 PISM Authors
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

#ifndef _PSCONSTANT_H_
#define _PSCONSTANT_H_

#include "PISMSurface.hh"
#include "iceModelVec.hh"
#include "PISMAtmosphere.hh"

//! \brief A class implementing a constant-in-time surface model.  Reads data
//! from a PISM input file.
/*!
This is model is just as simple as PSSimple, but it assumes results from a
surface processes model are already known.  But they are treated as constant in
time and they are read from the input file at the beginning of the PISM run.

Specifically, these two fields are read from the \c -i or \c -boot_file file:
  \li \c acab = ice-equivalent surface mass balance
  \li \c artm = ice fluid upper surface temperature.

This surface model does not use an atmosphere model at all, so the
\c attach_atmosphere_model() method is null.  Any choice of PISMAtmosphereModel
made using option \c -atmosphere is ignored.  This may be an advantage in coupler
code simplicity.

Note that a very minimal coupling of an existing atmosphere and surface processes
model to the ice dynamics core in PISM could be accomplished by using this
PSConstant class for relatively short ice dynamics runs, each of which starts by
reading the latest \c acab and \c artm fields supplied by the atmosphere and
surface processes model.
*/
class PSConstant : public PISMSurfaceModel {
public:
  PSConstant(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {};

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  { t = my_t; dt = my_dt; return 0; } // do nothing
  //! This surface model does not use an atmosphere model.
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  // Does not have an atmosphere model.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NetCDF3Wrapper &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual void add_vars_to_output(string keyword, set<string> &result);
protected:
  string input_file;
  IceModelVec2S acab, artm;
};

#endif /* _PSCONSTANT_H_ */
