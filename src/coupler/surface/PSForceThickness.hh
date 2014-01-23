// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#ifndef _PSFORCETHICKNESS_H_
#define _PSFORCETHICKNESS_H_

#include "PSModifier.hh"
#include "iceModelVec.hh"
#include "NCVariable.hh"

//! A class implementing a modified surface mass balance which forces
//! ice thickness to a given target by the end of the run.
class PSForceThickness : public PSModifier {
public:
  PSForceThickness(IceGrid &g, const PISMConfig &conf, PISMSurfaceModel *input);

  virtual ~PSForceThickness();
  virtual PetscErrorCode init(PISMVars &vars);
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input);
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict);
  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);
protected:
  std::string input_file;
  PetscReal alpha;
  IceModelVec2S *ice_thickness;	//!< current ice thickness produced by IceModel.
  IceModelVec2S target_thickness, ftt_mask;
  NCSpatialVariable climatic_mass_balance, climatic_mass_balance_original, ice_surface_temp;
private:
  PetscErrorCode allocate_PSForceThickness();
};

#endif /* _PSFORCETHICKNESS_H_ */
