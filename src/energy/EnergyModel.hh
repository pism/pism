/* Copyright (C) 2016, 2017, 2018, 2022, 2023 PISM Authors
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

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/CellType.hh"
#include <memory>

namespace pism {

namespace stressbalance {
class StressBalance;
}

namespace energy {

class Inputs {
public:
  Inputs();
  void check() const;

  const array::CellType *cell_type;
  const array::Scalar *basal_frictional_heating;
  const array::Scalar *basal_heat_flux;
  const array::Scalar1 *ice_thickness;
  const array::Scalar *surface_liquid_fraction;
  const array::Scalar *shelf_base_temp;
  const array::Scalar *surface_temp;
  const array::Scalar *till_water_thickness;

  const array::Array3D *volumetric_heating_rate;
  const array::Array3D *u3;
  const array::Array3D *v3;
  const array::Array3D *w3;

  // inputs used by regional models
  const array::Scalar *no_model_mask;
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
  EnergyModel(std::shared_ptr<const Grid> grid,
              std::shared_ptr<const stressbalance::StressBalance> stress_balance);

  void restart(const File &input_file, int record);

  /*! @brief Bootstrapping using heuristics. */
  /*!
   * Bootstrap by reading 2d fields (currently the basal melt rate) from a file and filling 3D
   * fields using heuristics.
   */
  void bootstrap(const File &input_file,
                 const array::Scalar &ice_thickness,
                 const array::Scalar &surface_temperature,
                 const array::Scalar &climatic_mass_balance,
                 const array::Scalar &basal_heat_flux);

  /*! @brief Initialize using formulas (for runs using synthetic data). */
  void initialize(const array::Scalar &basal_melt_rate,
                  const array::Scalar &ice_thickness,
                  const array::Scalar &surface_temperature,
                  const array::Scalar &climatic_mass_balance,
                  const array::Scalar &basal_heat_flux);

  void update(double t, double dt, const Inputs &inputs);

  const EnergyModelStats& stats() const;

  const array::Array3D & enthalpy() const;
  const array::Scalar & basal_melt_rate() const;

  const std::string& stdout_flags() const;
protected:

  virtual MaxTimestep max_timestep_impl(double t) const;

  virtual void restart_impl(const File &input_file, int record) = 0;

  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness,
                              const array::Scalar &surface_temperature,
                              const array::Scalar &climatic_mass_balance,
                              const array::Scalar &basal_heat_flux) = 0;

  virtual void initialize_impl(const array::Scalar &basal_melt_rate,
                               const array::Scalar &ice_thickness,
                               const array::Scalar &surface_temperature,
                               const array::Scalar &climatic_mass_balance,
                               const array::Scalar &basal_heat_flux) = 0;

  virtual void update_impl(double t, double dt, const Inputs &inputs) = 0;

  virtual void define_model_state_impl(const File &output) const = 0;
  virtual void write_model_state_impl(const File &output) const = 0;

  virtual DiagnosticList diagnostics_impl() const;
  virtual TSDiagnosticList ts_diagnostics_impl() const;

  /*! @brief Initialize enthalpy by reading it from a file, or by reading temperature and liquid
      water fraction, or by reading the temperature field alone. */
  void init_enthalpy(const File &input_file, bool regrid, int record);

  /*! @brief Regrid enthalpy from the -regrid_file. */
  void regrid_enthalpy();
protected:
  array::Array3D m_ice_enthalpy;
  array::Array3D m_work;
  array::Scalar m_basal_melt_rate;

  EnergyModelStats m_stats;

private:
  std::string m_stdout_flags;
  std::shared_ptr<const stressbalance::StressBalance> m_stress_balance;
};

/*!
 * Return true if the grid point (i,j) is near the margin of the ice.
 */
bool marginal(const array::Scalar1 &thickness, int i, int j, double threshold);

} // end of namespace energy
} // end of namespace pism


#endif /* ENERGYMODEL_H */
