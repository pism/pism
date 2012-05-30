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

#ifndef _PSDTFORCING_H_
#define _PSDTFORCING_H_

#include "PScalarForcing.hh"
#include "PISMSurface.hh"
#include "PSModifier.hh"

class PS_delta_T : public PScalarForcing<PISMSurfaceModel,PSModifier>
{
public:
  PS_delta_T(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel* in);
  virtual ~PS_delta_T() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
protected:
  NCSpatialVariable climatic_mass_balance, ice_surface_temp;
};

#endif /* _PSDTFORCING_H_ */
