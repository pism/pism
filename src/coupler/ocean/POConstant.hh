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

#ifndef _POCONSTANT_H_
#define _POCONSTANT_H_

#include "PISMOcean.hh"
#include "NCSpatialVariable.hh"

//! \brief A class implementing a constant (in terms of the ocean inputs) ocean
//! model. Uses configuration parameters for the sea level elevation and
//! sub-shelf heat flux.
class POConstant : public PISMOceanModel {
public:
  POConstant(IceGrid &g, const NCConfigVariable &conf);
  virtual ~POConstant() {}
  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt)
  { t = my_t; dt = my_dt; return 0; } // do nothing

  virtual PetscErrorCode sea_level_elevation(PetscReal &result);
  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
protected:
  IceModelVec2S *ice_thickness;	// is not owned by this class
  NCSpatialVariable shelfbmassflux, shelfbtemp;
  bool meltrate_set;
  PetscReal mymeltrate;
};

#endif /* _POCONSTANT_H_ */
