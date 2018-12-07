/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

#include <memory>               // shared_ptr

#include "pism/icemodel/IceModel.hh"

namespace pism {

namespace energy {
class CHSystem;
} // end of namespace energy

//! \brief A version of the PISM core class (IceModel) which knows about the
//! `no_model_mask` and its semantics.
class IceRegionalModel : public IceModel {
public:
  IceRegionalModel(IceGrid::Ptr g, Context::Ptr c);

  const energy::CHSystem* cryo_hydrologic_system() const;

protected:
  virtual void bootstrap_2d(const PIO &input_file);

  void allocate_geometry_evolution();
  void allocate_storage();
  void allocate_stressbalance();
  void allocate_basal_yield_stress();
  void allocate_energy_model();
  void model_state_setup();

  void energy_step();

  stressbalance::Inputs stress_balance_inputs();
  energy::Inputs energy_model_inputs();
  YieldStressInputs yield_stress_inputs();

  void init_diagnostics();

private:
  IceModelVec2Int m_no_model_mask;
  IceModelVec2S   m_usurf_stored;
  IceModelVec2S   m_thk_stored;

  std::shared_ptr<energy::CHSystem> m_ch_system;
  IceModelVec3::Ptr m_ch_warming_flux;
};

} // end of namespace pism

#endif /* _ICEREGIONALMODEL_H_ */
