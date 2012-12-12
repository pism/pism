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

#ifndef _PSSTUFFASANOMALY_H_
#define _PSSTUFFASANOMALY_H_

#include "PISMSurface.hh"
#include "PSModifier.hh"
#include "iceModelVec.hh"

//! \brief A surface modifier class applying its input as anomalies.
class PSStuffAsAnomaly : public PSModifier
{
public:
  PSStuffAsAnomaly(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel *input)
    : PSModifier(g, conf, input) {}
  virtual ~PSStuffAsAnomaly() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);

protected:
  IceModelVec2S mass_flux, mass_flux_0, mass_flux_input,
    temp, temp_0, temp_input;
};

#endif /* _PSSTUFFASANOMALY_H_ */
