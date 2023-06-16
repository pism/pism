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

#ifndef TEMPERATUREMODEL_VERIFICATION_H
#define TEMPERATUREMODEL_VERIFICATION_H

#include "pism/energy/TemperatureModel.hh"

namespace pism {
namespace energy {

/*! @brief Temperature-based energy balance model with verification-specific initialization code. */
class TemperatureModel_Verification : public TemperatureModel {
public:
  TemperatureModel_Verification(std::shared_ptr<const Grid> grid,
                                stressbalance::StressBalance *stress_balance,
                                int testname, bool bedrock_is_ice);

  void initialize_impl(const array::Scalar &basal_melt_rate,
                       const array::Scalar &ice_thickness,
                       const array::Scalar &surface_temperature,
                       const array::Scalar &climatic_mass_balance,
                       const array::Scalar &basal_heat_flux);
private:
  void initTestFG();
  void initTestsKO();

  int m_testname;
  bool m_bedrock_is_ice;
};

} // end of namespace energy
} // end of namespace pism

#endif /* TEMPERATUREMODEL_VERIFICATION_H */
