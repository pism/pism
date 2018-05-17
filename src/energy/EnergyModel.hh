/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#include "pism/util/Component.hh"

#include "pism/util/iceModelVec.hh"

namespace pism {

namespace stressbalance {
class StressBalance;
}

class IceModelVec2CellType;

namespace energy {

class Inputs {
public:
  Inputs();
  void check() const;

  const IceModelVec2CellType *cell_type;
  const IceModelVec2S *basal_frictional_heating;
  const IceModelVec2S *basal_heat_flux;
  const IceModelVec2S *ice_thickness;
  const IceModelVec2S *surface_liquid_fraction;
  const IceModelVec2S *shelf_base_temp;
  const IceModelVec2S *surface_temp;
  const IceModelVec2S *till_water_thickness;

  const IceModelVec3 *volumetric_heating_rate;
  const IceModelVec3 *u3;
  const IceModelVec3 *v3;
  const IceModelVec3 *w3;

  // inputs used by regional models
  const IceModelVec2Int *no_model_mask;
};

class EnergyModelStats {
public:
  EnergyModelStats();

  EnergyModelStats& operator+=(const EnergyModelStats &other);

  void sum(MPI_Comm com);

  unsigned int bulge_counter;
  unsigned int reduced_accuracy_counter;
  unsigned int low_temperature_counter;
  double liquified_ice_volume;
};

class EnergyModel : public Component {
public:
  EnergyModel(IceGrid::ConstPtr grid, stressbalance::StressBalance *stress_balance);

  void restart(const PIO &input_file, int record);

  /*! @brief Bootstrapping using heuristics. */
  /*!
   * Bootstrap by reading 2d fields (currently the basal melt rate) from a file and filling 3D
   * fields using heuristics.
   */
  void bootstrap(const PIO &input_file,
                 const IceModelVec2S &ice_thickness,
                 const IceModelVec2S &surface_temperature,
                 const IceModelVec2S &climatic_mass_balance,
                 const IceModelVec2S &basal_heat_flux);

  /*! @brief Initialize using formulas (for runs using synthetic data). */
  void initialize(const IceModelVec2S &basal_melt_rate,
                  const IceModelVec2S &ice_thickness,
                  const IceModelVec2S &surface_temperature,
                  const IceModelVec2S &climatic_mass_balance,
                  const IceModelVec2S &basal_heat_flux);

  void update(double t, double dt, const Inputs &inputs);

  const EnergyModelStats& stats() const;

  const IceModelVec3 & enthalpy() const;
  const IceModelVec2S & basal_melt_rate() const;

  const std::string& stdout_flags() const;
protected:

  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void restart_impl(const PIO &input_file, int record) = 0;

  virtual void bootstrap_impl(const PIO &input_file,
                              const IceModelVec2S &ice_thickness,
                              const IceModelVec2S &surface_temperature,
                              const IceModelVec2S &climatic_mass_balance,
                              const IceModelVec2S &basal_heat_flux) = 0;

  virtual void initialize_impl(const IceModelVec2S &basal_melt_rate,
                               const IceModelVec2S &ice_thickness,
                               const IceModelVec2S &surface_temperature,
                               const IceModelVec2S &climatic_mass_balance,
                               const IceModelVec2S &basal_heat_flux) = 0;

  virtual void update_impl(double t, double dt, const Inputs &inputs) = 0;

  virtual void define_model_state_impl(const PIO &output) const = 0;
  virtual void write_model_state_impl(const PIO &output) const = 0;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  /*! @brief Initialize enthalpy by reading it from a file, or by reading temperature and liquid
      water fraction, or by reading the temperature field alone. */
  void init_enthalpy(const PIO &input_file, bool regrid, int record);

  /*! @brief Regrid enthalpy from the -regrid_file. */
  void regrid_enthalpy();
protected:
  IceModelVec3 m_ice_enthalpy;
  IceModelVec3 m_work;
  IceModelVec2S m_basal_melt_rate;

  EnergyModelStats m_stats;

private:
  std::string m_stdout_flags;
  stressbalance::StressBalance *m_stress_balance;
};

/*!
 * Return true if the grid point (i,j) is near the margin of the ice.
 */
bool marginal(const IceModelVec2S &thickness, int i, int j, double threshold);

} // end of namespace energy
} // end of namespace pism


#endif /* ENERGYMODEL_H */
