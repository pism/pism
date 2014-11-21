// Copyright (C) 2011, 2014 PISM Authors
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

#ifndef _PSMODIFIER_H_
#define _PSMODIFIER_H_

#include "PISMSurface.hh"

namespace pism {

//! \brief A base class for mechanisms which modify the results of a surface
//! processes model (an instance of SurfaceModel) before they reach the ice.
/*! 
  Frequently ice sheet models are driven by a "basic" surface model plus "forcings".
  This modifier class allows the implementations of forcings which alter the 
  results of the surface processes model.  That is, if the atmospheric inputs 
  are already dealt-with, and a basic surface processes model is in use which 
  generates surface mass balance and ice upper surface temperature, then instances
  of this PSModifier class can be used to modify the surface mass balance and ice
  upper surface temperature "just before" it gets to the ice itself.
*/
class PSModifier : public Modifier<SurfaceModel>
{
public:
  PSModifier(IceGrid &g, const Config &conf, SurfaceModel* in)
    : Modifier<SurfaceModel>(g, conf, in) {}
  virtual ~PSModifier() {}

  virtual void attach_atmosphere_model(AtmosphereModel *in) {
    if (input_model != NULL) {
      input_model->attach_atmosphere_model(in);
    }
  }

  virtual void ice_surface_mass_flux(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->ice_surface_mass_flux(result);
    }
  }

  virtual void ice_surface_temperature(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->ice_surface_temperature(result);
    }
  }

  virtual void ice_surface_liquid_water_fraction(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->ice_surface_liquid_water_fraction(result);
    }
  }

  virtual void mass_held_in_surface_layer(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->mass_held_in_surface_layer(result);
    }
  }

  virtual void surface_layer_thickness(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      input_model->surface_layer_thickness(result);
    }
  }
};

} // end of namespace pism

#endif /* _PSMODIFIER_H_ */
