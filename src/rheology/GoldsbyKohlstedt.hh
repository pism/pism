/* Copyright (C) 2015, 2016, 2021 PISM Authors
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

#ifndef _GOLDSBYKOHLSTEDT_H_
#define _GOLDSBYKOHLSTEDT_H_

#include "FlowLaw.hh"

namespace pism {
namespace rheology {

// Hybrid (Goldsby-Kohlstedt/Glen) ice flow law

struct GKparts {
  double eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};

//! A hybrid of Goldsby-Kohlstedt (2001) ice (constitutive form) and Paterson-Budd (1982)-Glen (viscosity form).
/*!
  Each FlowLaw has both a forward flow law in "constitutive law" form ("flow_from_temp()") and an
  inverted-and-vertically-integrated flow law ("effective_viscosity()").  Only the
  former form of the flow law is known for Goldsby-Kohlstedt.  If one can
  invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
  derived class.
*/
class GoldsbyKohlstedt : public FlowLaw {
public:
  GoldsbyKohlstedt(const std::string &prefix,
                   const Config &config,
                   EnthalpyConverter::Ptr EC);
protected:
  virtual double flow_impl(double stress, double E,
                           double pressure, double grainsize) const;

  // NB! not virtual
  double softness_impl(double E, double p) const __attribute__((noreturn));
  double hardness_impl(double E, double p) const;
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;
  GKparts flowParts(double stress, double temp, double pressure) const;

  double  m_V_act_vol,  m_d_grain_size,
  //--- diffusional flow ---
    m_diff_crit_temp, m_diff_V_m, m_diff_D_0v, m_diff_Q_v, m_diff_D_0b, m_diff_Q_b, m_diff_delta,
  //--- dislocation creep ---
    m_disl_crit_temp, m_disl_A_cold, m_disl_A_warm, m_disl_n, m_disl_Q_cold, m_disl_Q_warm,
  //--- easy slip (basal) ---
    m_basal_A, m_basal_n, m_basal_Q,
  //--- grain boundary sliding ---
    m_gbs_crit_temp, m_gbs_A_cold, m_gbs_A_warm, m_gbs_n, m_gbs_Q_cold,
    m_p_grain_sz_exp, m_gbs_Q_warm;
};

//! Derived class of GoldsbyKohlstedt for testing purposes only.
/*!
  GoldsbyKohlstedtStripped is a simplification of Goldsby-Kohlstedt. Compare to that
  used in Peltier et al 2000, which is even simpler.
*/
class GoldsbyKohlstedtStripped : public GoldsbyKohlstedt {
public:
  GoldsbyKohlstedtStripped(const std::string &prefix,
                           const Config &config, EnthalpyConverter::Ptr EC);

protected:
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;

  double m_d_grain_size_stripped;
};

} // end of namespace rheology
} // end of namespace pism


#endif /* _GOLDSBYKOHLSTEDT_H_ */
