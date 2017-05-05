/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

#ifndef _ICEREGIONALMODEL_H_
#define _ICEREGIONALMODEL_H_

#include "base/iceModel.hh"

namespace pism {

//! \brief A version of the PISM core class (IceModel) which knows about the
//! `no_model_mask` and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid::Ptr g, Context::Ptr c);
protected:
  virtual void bootstrap_2d(const PIO &input_file);

  virtual void allocate_geometry_evolution();
  virtual void allocate_storage();
  virtual void allocate_stressbalance();
  virtual void allocate_basal_yield_stress();
  virtual void allocate_energy_model();
private:
  IceModelVec2Int m_no_model_mask;
  IceModelVec2S   m_usurf_stored;
  IceModelVec2S   m_thk_stored;
};

} // end of namespace pism

#endif /* _ICEREGIONALMODEL_H_ */
