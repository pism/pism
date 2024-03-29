/* Copyright (C) 2016, 2017, 2023 PISM Authors
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

#ifndef ENTHALPYMODEL_REGIONAL_H
#define ENTHALPYMODEL_REGIONAL_H

#include "pism/energy/EnthalpyModel.hh"

namespace pism {
namespace energy {

/*! @brief The enthalpy-based energy balance model for regional runs. */
class EnthalpyModel_Regional : public EnthalpyModel {
public:
  EnthalpyModel_Regional(std::shared_ptr<const Grid> grid,
                         std::shared_ptr<const stressbalance::StressBalance> stress_balance);

protected:
  virtual void restart_impl(const File &input_file, int record);

  virtual void bootstrap_impl(const File &input_file,
                              const array::Scalar &ice_thickness,
                              const array::Scalar &surface_temperature,
                              const array::Scalar &climatic_mass_balance,
                              const array::Scalar &basal_heat_flux);

  virtual void initialize_impl(const array::Scalar &basal_melt_rate,
                               const array::Scalar &ice_thickness,
                               const array::Scalar &surface_temperature,
                               const array::Scalar &climatic_mass_balance,
                               const array::Scalar &basal_heat_flux);

  void update_impl(double t, double dt, const Inputs &inputs);

  array::Scalar *m_no_model_mask;
  array::Scalar m_basal_melt_rate_stored;
};

} // end of namespace energy
} // end of namespace pism


#endif /* ENTHALPYMODEL_REGIONAL_H */
