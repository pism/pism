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

#ifndef _PS_DELTA_T_H_
#define _PS_DELTA_T_H_

#include "PScalarForcing.hh"
#include "PISMSurface.hh"
#include "PSModifier.hh"

namespace pism {

/** @brief Implements the scalar temperature offsets for the ice
 * surface temperature.
 *
 * Other fields are passed through without change.
 */
class PS_delta_T : public PScalarForcing<SurfaceModel,PSModifier>
{
public:
  PS_delta_T(IceGrid &g, const Config &conf, SurfaceModel* in);
  virtual ~PS_delta_T();

  virtual void init(Vars &vars);

  virtual void ice_surface_temperature(IceModelVec2S &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc,
                                          IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output(const std::string &keyword, std::set<std::string> &result);
protected:
  NCSpatialVariable climatic_mass_balance, //!< climatic mass balance attributes
    ice_surface_temp;                      //!< ice surface temperature attributes
private:
  //! Allocate internal objects. Called from the constructor.
  PetscErrorCode allocate_PS_delta_T();
};

} // end of namespace pism

#endif /* _PS_DELTA_T_H_ */
