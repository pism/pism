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

namespace ocean {

struct TocBox1 {
  bool failed;
  double value;
};

class PicoPhysics {
public:
  PicoPhysics(const Config &config);

  double pressure(double ice_thickness) const;
  double T_star(double salinity, double temperature, double pressure) const;

  TocBox1 Toc_box1(double area, double T_star, double Soc_box0, double Toc_box0) const;
  double Soc_box1(double Toc_box0, double Soc_box0, double Toc) const;

  double Toc(double box_area, double temperature, double T_star, double overturning, double salinity) const;

  double Soc(double salinity, double temperature, double Toc) const;

  double theta_pm(double salinity, double pressure) const;
  double T_pm(double salinity, double pressure) const;

  double melt_rate(double pm_point, double Toc) const;

  double melt_rate_beckmann_goosse(double pot_pm_point, double Toc) const;

  double overturning(double Soc_box0, double Soc, double Toc_box0, double Toc) const;

  double gamma_T() const;
  double overturning_coeff() const;
  double T_dummy() const;
  double S_dummy() const;
  double ice_density() const;
  double continental_shelf_depth() const;

private:
  double p_coeff(double g1, double s1) const;

  double m_gamma_T, m_overturning_coeff, m_T_dummy, m_S_dummy;
  double m_ice_density, m_continental_shelf_depth;

  double m_earth_grav, m_sea_water_density, m_rho_star, m_nu, m_latentHeat, m_c_p_ocean, m_alpha, m_beta;

  double m_lambda;

  // coefficients of the parameterization of the potential temperature
  double m_a_pot, m_b_pot, m_c_pot;

  // coefficients of the parameterization of the in situ temperature
  double m_a_in_situ, m_b_in_situ, m_c_in_situ;

  double m_meltFactor;
};

} // end of namespace ocean
} // end of namespace pism
