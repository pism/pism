/* Copyright (C) 2016 PISM Authors
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

#ifndef AGEMODEL_H
#define AGEMODEL_H

#include "base/util/IceModelVec.hh"
#include "base/util/PISMComponent.hh"
#include "base/stressbalance/PISMStressBalance.hh"

namespace pism {

class AgeModelInputs {
public:
  AgeModelInputs();
  void check() const;

  const IceModelVec2S *ice_thickness;
  const IceModelVec3 *u3;
  const IceModelVec3 *v3;
  const IceModelVec3 *w3;
};

class AgeModel : public Component_TS {
public:
  AgeModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

  using Component_TS::update;
  void update(double t, double dt, const AgeModelInputs &inputs);

  void init(const InputOptions &opts);

protected:
  void init_impl(const InputOptions &opts);
  MaxTimestep max_timestep_impl(double t);
  void update_impl(double t, double dt);

  IceModelVec3 m_ice_age;
  IceModelVec3 m_work;
  stressbalance::StressBalance *m_stress_balance;
};

} // end of namespace pism


#endif /* AGEMODEL_H */
