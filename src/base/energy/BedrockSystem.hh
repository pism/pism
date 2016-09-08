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

#ifndef BEDROCKSYSTEM_H
#define BEDROCKSYSTEM_H

namespace pism {
namespace energy {

class BedrockSystem : public columnSystemCtx {
public:
  BedrockSystem(const std::string &prefix,
                double dt, double dz,
                const Config &config,
                const IceModelVec3Custom &temp);
  ~BedrockSystem();

  void init(int i, int j);

  void set_top_surface_dirichlet(double T);
  void set_bottom_surface_neumann(double dT);
  void set_bottom_surface_heat_flux(double Q);

  void solve(std::vector<double> &result);

protected:
  double m_D_bottom, m_U_bottom, m_B_bottom;
  double m_D_top, m_U_top, m_B_top;
};

} // end of namespace energy
} // end of namespace pism

#endif /* BEDROCKSYSTEM_H */
