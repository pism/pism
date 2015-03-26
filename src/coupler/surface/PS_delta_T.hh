// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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
namespace surface {

/** @brief Implements the scalar temperature offsets for the ice
 * surface temperature.
 *
 * Other fields are passed through without change.
 */
class Delta_T : public PScalarForcing<SurfaceModel,SurfaceModifier>
{
public:
  Delta_T(const IceGrid &g, SurfaceModel* in);
  virtual ~Delta_T();

  virtual void init();

  virtual void ice_surface_temperature(IceModelVec2S &result);

protected:
  virtual MaxTimestep max_timestep_impl(double t);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype);
protected:
  SpatialVariableMetadata climatic_mass_balance, //!< climatic mass balance attributes
    ice_surface_temp;                      //!< ice surface temperature attributes
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PS_DELTA_T_H_ */
