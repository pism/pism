// Copyright (C) 2011, 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _POMODIFIER_H_
#define _POMODIFIER_H_

#include "coupler/PISMOcean.hh"
#include "base/util/Modifier.hh"

namespace pism {
namespace ocean {
class OceanModifier : public Modifier<OceanModel>
{
public:
  OceanModifier(IceGrid::ConstPtr g, OceanModel* in)
    : Modifier<OceanModel>(g, in) {}
  virtual ~OceanModifier() {}

protected:
  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result)
  {
    m_input_model->melange_back_pressure_fraction(result);
  }
  virtual void shelf_base_temperature_impl(IceModelVec2S &result)
  {
    m_input_model->shelf_base_temperature(result);
  }

  virtual void sea_level_elevation_impl(double &result)
  {
    result = m_input_model->sea_level_elevation();
  }
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result)
  {
    m_input_model->shelf_base_mass_flux(result);
  }
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POMODIFIER_H_ */
