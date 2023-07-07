// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PO_RUNOFF_SMB_H_
#define _PO_RUNOFF_SMB_H_

#include "pism/coupler/OceanModel.hh"

namespace pism {

class ScalarForcing;

namespace ocean {

/** Modify the shelf base mass flux using a function of air temperature changes.
 *
 * This modifier multiplies the input shelf base mass flux by a factor
 * @f$ F @f$, which depends on the change in the global average
 * near-surface air temperature =delta_T= read from an input file.
 *
 * The parameterization described below is inspired by [@ref Xu2013]
 * with the relationship between air temperatures and subglacial
 * runoff flux parameterized using a linear function (a fit to
 * results of GCM simulations).
 *
 * @f[ F(\Delta T) = 1 + B \times (a \Delta T)^\alpha \times (\Delta T)^\beta. @f]
 *
 * Here @f$ a @f$, @f$ B @f$, @f$ \alpha @f$, and @f$ \beta @f$ are constants.
 *
 * The paper [@ref Xu2013] approximates the sub-shelf melt rate as a
 * function of *ocean* temperature above the freezing point and the
 * subglacial runoff.
 *
 * We assume that the lag between the air and ocean temperatures is
 * negligible and a change in air temperature is directly translated
 * into a change in ocean temperature.
 */
class Runoff_SMB : public OceanModel
{
public:
  Runoff_SMB(std::shared_ptr<const Grid> g, std::shared_ptr<OceanModel> in);
  virtual ~Runoff_SMB();

private:
  void init_impl(const Geometry &geometry);

  void update_impl(const Geometry &geometry, double t, double dt);

  void mass_flux(double delta_T, array::Scalar &result) const;

  // @brief constant in the parameterization of the subglacial
  // runoff flux as a function of air temperature
  double m_temp_to_runoff_a;

  // Constants in the parameterization of the sub-shelf melt rate.
  double m_runoff_to_ocean_melt_b;
  double m_runoff_to_ocean_melt_power_alpha;
  double m_runoff_to_ocean_melt_power_beta;

  std::shared_ptr<array::Scalar> m_shelf_base_mass_flux;
  std::unique_ptr<ScalarForcing> m_forcing;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _PO_RUNOFF_SMB_H_ */
