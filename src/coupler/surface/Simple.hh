// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PSSIMPLE_H_
#define _PSSIMPLE_H_

#include "pism/coupler/SurfaceModel.hh"

namespace pism {
namespace surface {

//! \brief A class implementing a primitive surface model.
/*! 
  This is an "invisible" surface processes model which "passes through"
  information from the atmosphere above directly to the ice below the surface
  layers.  It implements two modeling choices:
  \li accumulation which is obtained from an atmosphere model is interpreted
  as surface mass flux;
  \li mean-annual near-surface air temperature is interpreted as instantaneous
  temperature of the ice at the ice surface.

  The second choice means that the upper boundary condition of the conservation of
  energy scheme for the ice fluid is exactly the 2m air temperature.
*/
class Simple : public SurfaceModel {
public:
  Simple(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> atmosphere);
protected:
  virtual void init_impl(const Geometry &geometry);
  virtual void update_impl(const Geometry &geometry, double t, double dt);

  virtual const IceModelVec2S& mass_flux_impl() const;
  virtual const IceModelVec2S& temperature_impl() const;

  IceModelVec2S::Ptr m_mass_flux;
  IceModelVec2S::Ptr m_temperature;
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSSIMPLE_H_ */
