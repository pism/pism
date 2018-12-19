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

namespace pism {

class Config;

namespace frontalmelt {

class FrontalMeltPhysics {
public:
  FrontalMeltPhysics(const Config &config);

  double frontal_melt_from_undercutting(double ice_thickness,
                                        double discharge_flux,
                                        double potential_temperature) const;

  double frontal_melt_from_ismip6(double ice_thickness,
                                        double discharge_flux,
                                        double potential_temperature) const;

private:
  double m_A, m_B, m_alpha, m_beta;
};

} // end of namespace frontalmelt
} // end of namespace pism
