/* Copyright (C) 2018 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef COMPLETEOCEANMODEL_H
#define COMPLETEOCEANMODEL_H

#include "pism/coupler/OceanModel.hh"

namespace pism {
namespace ocean {

/*!
 * Base lass for ocean models with dedicated storage for output fields.
 *
 * All ocean models have storage for melange back pressure. (All but one set it to zero.)
 */
class CompleteOceanModel : public OceanModel {
public:
  // "modifier" constructor
  CompleteOceanModel(IceGrid::ConstPtr g, std::shared_ptr<OceanModel> input);
  // "model" constructor
  CompleteOceanModel(IceGrid::ConstPtr g);

  virtual ~CompleteOceanModel();
protected:
  virtual const IceModelVec2S& shelf_base_temperature_impl() const;
  virtual const IceModelVec2S& shelf_base_mass_flux_impl() const;
  // getter for melange_back_pressure_fraction is inherited from OceanModel

  // storage for melange_back_pressure_fraction is inherited from OceanModel
  IceModelVec2S::Ptr m_shelf_base_temperature;
  IceModelVec2S::Ptr m_shelf_base_mass_flux;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* COMPLETEOCEANMODEL_H */
