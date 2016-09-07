// Copyright (C) 2011, 2014, 2015, 2016 PISM Authors
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

#include "coupler/PISMSurface.hh"
#include "base/util/Modifier.hh"

namespace pism {
namespace surface {

//! \brief A base class for mechanisms which modify the results of a surface
//! processes model (an instance of SurfaceModel) before they reach the ice.
/*! 
  Frequently ice sheet models are driven by a "basic" surface model plus "forcings".
  This modifier class allows the implementations of forcings which alter the 
  results of the surface processes model.  That is, if the atmospheric inputs 
  are already dealt-with, and a basic surface processes model is in use which 
  generates surface mass balance and ice upper surface temperature, then instances
  of this SurfaceModifier class can be used to modify the surface mass balance and ice
  upper surface temperature "just before" it gets to the ice itself.
*/
class SurfaceModifier : public Modifier<SurfaceModel>
{
public:
  SurfaceModifier(IceGrid::ConstPtr g, SurfaceModel* in)
    : Modifier<SurfaceModel>(g, in) {}
  virtual ~SurfaceModifier() {}
protected:
  virtual void attach_atmosphere_model_impl(atmosphere::AtmosphereModel *in) {
    if (m_input_model != NULL) {
      m_input_model->attach_atmosphere_model(in);
    }
  }

  virtual void ice_surface_mass_flux_impl(IceModelVec2S &result)
  {
    if (m_input_model != NULL) {
      m_input_model->ice_surface_mass_flux(result);
    }
  }

  virtual void ice_surface_temperature_impl(IceModelVec2S &result)
  {
    if (m_input_model != NULL) {
      m_input_model->ice_surface_temperature(result);
    }
  }

  virtual void ice_surface_liquid_water_fraction_impl(IceModelVec2S &result)
  {
    if (m_input_model != NULL) {
      m_input_model->ice_surface_liquid_water_fraction(result);
    }
  }

  virtual void mass_held_in_surface_layer_impl(IceModelVec2S &result)
  {
    if (m_input_model != NULL) {
      m_input_model->mass_held_in_surface_layer(result);
    }
  }

  virtual void surface_layer_thickness_impl(IceModelVec2S &result)
  {
    if (m_input_model != NULL) {
      m_input_model->surface_layer_thickness(result);
    }
  }
};

} // end of namespace surface
} // end of namespace pism

#endif /* _PSMODIFIER_H_ */
