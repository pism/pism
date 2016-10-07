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

#ifndef ENERGYMODEL_H
#define ENERGYMODEL_H

#include "base/util/PISMComponent.hh"

#include "base/util/iceModelVec.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {

namespace stressbalance {
class StressBalance;
}

namespace energy {

class EnergyModelInputs {
public:
  EnergyModelInputs();
  void check() const;

  const IceModelVec2CellType *cell_type;
  const IceModelVec2S *basal_frictional_heating;
  const IceModelVec2S *basal_heat_flux;
  const IceModelVec2S *ice_thickness;
  const IceModelVec2S *surface_liquid_fraction;
  const IceModelVec2S *shelf_base_temp;
  const IceModelVec2S *surface_temp;
  const IceModelVec2S *till_water_thickness;

  const IceModelVec3 *strain_heating3;
  const IceModelVec3 *u3;
  const IceModelVec3 *v3;
  const IceModelVec3 *w3;
};

class EnergyModelStats {
public:
  EnergyModelStats();

  unsigned int bulge_counter;
  unsigned int reduced_accuracy_counter;
  unsigned int low_temperature_counter;
  double liquified_ice_volume;
};

class EnergyModel : public Component_TS {
public:
  EnergyModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

  void init(const InputOptions &opts);

  using Component_TS::update;
  void update(double t, double dt, const EnergyModelInputs &inputs);

  const EnergyModelStats& stats() const;
  const IceModelVec3 & enthalpy() const;
  const IceModelVec2S & basal_melt_rate() const;
protected:
  void update_impl(double t, double dt);
  MaxTimestep max_timestep_impl(double t) const;

  virtual void init_impl(const InputOptions &opts) = 0;
  virtual void update_impl(double t, double dt, const EnergyModelInputs &inputs) = 0;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  void init_enthalpy(const PIO &input_file, bool regrid, int record);
protected:
  IceModelVec3 m_ice_enthalpy;
  IceModelVec3 m_work;
  IceModelVec2S m_basal_melt_rate;

  EnergyModelStats m_stats;

  stressbalance::StressBalance *m_stress_balance;
};

class EnthalpyModel : public EnergyModel {
public:
  EnthalpyModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

protected:
  void init_impl(const InputOptions &opts);
  void update_impl(double t, double dt, const EnergyModelInputs &inputs);
};

class TemperatureModel : public EnergyModel {
public:
  TemperatureModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

  const IceModelVec3 & temperature() const;
protected:
  void init_impl(const InputOptions &opts);
  void update_impl(double t, double dt, const EnergyModelInputs &inputs);

  void define_model_state_impl(const PIO &output) const;
  void write_model_state_impl(const PIO &output) const;

  void column_drainage(const double rho, const double c, const double L,
                       const double z, const double dz,
                       double *Texcess, double *bwat) const;

  IceModelVec3 m_ice_temperature;
};

} // end of namespace energy
} // end of namespace pism


#endif /* ENERGYMODEL_H */
