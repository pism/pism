/* Copyright (C) 2016, 2017 PISM Authors
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

#include "pism/util/iceModelVec.hh"
#include "pism/util/Component.hh"
#include "pism/stressbalance/StressBalance.hh"

namespace pism {

class AgeModelInputs {
public:
  AgeModelInputs();
  AgeModelInputs(const IceModelVec2S *ice_thickness,
                 const IceModelVec3 *u3,
                 const IceModelVec3 *v3,
                 const IceModelVec3 *w3);
  void check() const;

  const IceModelVec2S *ice_thickness;
  const IceModelVec3 *u3;
  const IceModelVec3 *v3;
  const IceModelVec3 *w3;
};

class AgeModel : public Component {
public:
  AgeModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

  void update(double t, double dt, const AgeModelInputs &inputs);

  void init(const InputOptions &opts);

  const IceModelVec3 & age() const;
protected:
  MaxTimestep max_timestep_impl(double t) const;
  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  IceModelVec3 m_ice_age;
  IceModelVec3 m_work;
  stressbalance::StressBalance *m_stress_balance;
};

} // end of namespace pism


#endif /* AGEMODEL_H */
