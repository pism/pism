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

#ifndef CHSYSTEM_H
#define CHSYSTEM_H

#include "EnergyModel.hh"

namespace pism {
namespace energy {

class CHSystem : public EnergyModel {
public:
  CHSystem(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);
  virtual ~CHSystem();

protected:
  void restart_impl(const PIO &input_file, int record);

  void bootstrap_impl(const PIO &input_file,
                      const IceModelVec2S &ice_thickness,
                      const IceModelVec2S &surface_temperature,
                      const IceModelVec2S &climatic_mass_balance,
                      const IceModelVec2S &basal_heat_flux);

  void initialize_impl(const IceModelVec2S &basal_melt_rate,
                       const IceModelVec2S &ice_thickness,
                       const IceModelVec2S &surface_temperature,
                       const IceModelVec2S &climatic_mass_balance,
                       const IceModelVec2S &basal_heat_flux);

  void update_impl(double t, double dt, const Inputs &inputs);

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;
};

void cryo_hydrologic_warming_flux(double k,
                                  double R,
                                  const IceModelVec2S &ice_thickness,
                                  const IceModelVec3 &ice_enthalpy,
                                  const IceModelVec3 &ch_enthalpy,
                                  IceModelVec3 &result);

} // end of namespace energy
} // end of namespace pism

#endif /* CHSYSTEM_H */
