// Copyright (C) 2011, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "pism/coupler/OceanModel.hh"
#include "pism/util/Modifier.hh"

namespace pism {
namespace ocean {
class OceanModifier : public Modifier<OceanModel>
{
public:
  OceanModifier(IceGrid::ConstPtr g, OceanModel* in)
    : Modifier<OceanModel>(g, in) {}
  virtual ~OceanModifier() {}

protected:
  void update_impl(double t, double dt) {
    m_input_model->update(t, dt);

    m_sea_level = m_input_model->sea_level_elevation();
    m_melange_back_pressure_fraction.copy_from(m_input_model->melange_back_pressure_fraction());
    m_shelf_base_temperature.copy_from(m_input_model->shelf_base_temperature());
    m_shelf_base_mass_flux.copy_from(m_input_model->shelf_base_mass_flux());
  }
};

} // end of namespace ocean
} // end of namespace pism
#endif /* _POMODIFIER_H_ */
