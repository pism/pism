/* Copyright (C) 2015, 2016 PISM Authors
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
  virtual void restart_2d(const PIO &input_file, unsigned int record);
  virtual void model_state_setup();
  virtual void createVecs();
  virtual void allocate_stressbalance();
  virtual void allocate_basal_yield_stress();
  virtual void massContExplicitStep(double dt,
                                    const IceModelVec2Stag &diffusive_flux,
                                    const IceModelVec2V &advective_velocity);
  virtual void cell_interface_fluxes(int i, int j,
                                     StarStencil<int> cell_type,
                                     StarStencil<int> bc_mask,
                                     StarStencil<Vector2> bc_velocity,
                                     StarStencil<Vector2> in_SSA_velocity,
                                     StarStencil<double> in_SIA_flux,
                                     StarStencil<double> &out_SSA_velocity,
                                     StarStencil<double> &out_SIA_flux);
private:
  IceModelVec2Int m_no_model_mask;
  IceModelVec2S   m_usurf_stored, m_thk_stored;
};

} // end of namespace pism

#endif /* _ICEREGIONALMODEL_H_ */
