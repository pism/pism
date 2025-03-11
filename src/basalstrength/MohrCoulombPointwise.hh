/* Copyright (C) 2019, 2025 PISM Authors
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
#ifndef MOHRCOULOMBPOINTWISE_H
#define MOHRCOULOMBPOINTWISE_H

#include "pism/util/Config.hh"

namespace pism {

/*!
 * Implementation of the Mohr-Coulomb basal yield stress model *at a point*.
 */
class MohrCoulombPointwise {
public:
  MohrCoulombPointwise(std::shared_ptr<const Config> config);

  /*!
   * Compute basal yield stress
   *
   * @param[in] delta fraction of overburden pressure
   * @param[in] P_overburden overburden pressure (Pa)
   * @param[in] water_thickness till water thickness (m)
   * @param[in] phi till friction angle (degrees)
   *
   * returns basal yield stress in Pascal
   */
  double yield_stress(double delta,
                      double P_overburden,
                      double water_thickness,
                      double phi) const;

  /*!
   * Inverse of `yield_stress()`.
   *
   * @param[in] delta fraction of overburden pressure
   * @param[in] P_overburden overburden pressure (Pa)
   * @param[in] water_thickness till water thickness
   * @param[in] yield_stress basal yield stress (Pa)
   *
   * returns till friction angle in degrees
   */
  double till_friction_angle(double delta,
                             double P_overburden,
                             double water_thickness,
                             double yield_stress) const;

  /*!
   * Compute effective pressure on till
   *
   * Used in `yield_stress()` and `till_friction_angle()`.
   *
   * @param[in] delta fraction of overburden pressure
   * @param[in] P_overburden overburden pressure (Pa)
   * @param[in] water_thickness till water thickness
   */
  double effective_pressure(double delta,
                            double P_overburden,
                            double water_thickness) const;
private:
  //! Maximum till water thickness
  double m_W_till_max;
  //! Cohesion of till
  double m_till_cohesion;
  //! Reference effective pressure
  double m_reference_effective_pressure;
  //! Reference void ratio
  double m_reference_void_ratio;
  //! Coefficient of compressibility of till
  double m_compressibility_coefficient;
};

} // end of namespace pism


#endif /* MOHRCOULOMBPOINTWISE_H */
