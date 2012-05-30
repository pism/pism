// Copyright (C) 2011, 2012 Andy Aschwanden and Constantine Khroulev
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

#ifndef _PSELEVATION_H_
#define _PSELEVATION_H_

#include "PISMSurface.hh"
#include "iceModelVec2T.hh"
#include "PISMAtmosphere.hh"

//! \brief A class implementing a elevation-dependent temperature and mass balance model.

class PSElevation : public PISMSurfaceModel {
public:
  PSElevation(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {};

  virtual PetscErrorCode init(PISMVars &vars);
  //! This surface model does not use an atmosphere model.
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  // Does not have an atmosphere model.
  virtual void get_diagnostics(map<string, PISMDiagnostic*> &/*dict*/) {}

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  { t = my_t; dt = my_dt; return 0; } // do nothing
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
protected:
  NCSpatialVariable climatic_mass_balance, ice_surface_temp;
  IceModelVec2S *usurf;
  PetscReal T_min, T_max, z_T_min, z_T_max,
    m_min, m_max, m_limit_min, m_limit_max,
    z_m_min, z_ELA, z_m_max;
};

#endif /* _PSELEVATION_H_ */
