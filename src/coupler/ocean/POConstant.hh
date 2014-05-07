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

#ifndef _POCONSTANT_H_
#define _POCONSTANT_H_

#include "PISMOcean.hh"
#include "NCVariable.hh"

namespace pism {

//! \brief A class implementing a constant (in terms of the ocean inputs) ocean
//! model. Uses configuration parameters for the sea level elevation and
//! sub-shelf heat flux.
class POConstant : public OceanModel {
public:
  POConstant(IceGrid &g, const Config &conf);
  virtual ~POConstant() {}

  virtual PetscErrorCode init(Vars &vars);

  virtual PetscErrorCode update(double my_t, double my_dt)
  { m_t = my_t; m_dt = my_dt; return 0; } // do nothing

  virtual PetscErrorCode sea_level_elevation(double &result);
  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
  virtual PetscErrorCode define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(const std::set<std::string> &vars, const PIO &nc);
protected:
  IceModelVec2S *ice_thickness; // is not owned by this class
  NCSpatialVariable shelfbmassflux, shelfbtemp;
  bool meltrate_set;
  double mymeltrate;
private:
  PetscErrorCode allocate_POConstant();
};

} // end of namespace pism

#endif /* _POCONSTANT_H_ */
