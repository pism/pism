// Copyright (C) 2004-2015 Jed Brown, Ed Bueler, and Constantine Khroulev
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

#include "flowlaws.hh"
#include "base/util/pism_const.hh"
#include "base/enthalpyConverter.hh"
#include "base/util/pism_options.hh"
#include "base/util/iceModelVec.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"

#include "base/util/error_handling.hh"

namespace pism {
namespace rheology {

bool FlowLawUsesGrainSize(FlowLaw *flow_law) {
  static const double gs[] = {1e-4, 1e-3, 1e-2, 1}, s=1e4, E=400000, p=1e6;
  double ref = flow_law->flow(s, E, p, gs[0]);
  for (int i=1; i<4; i++) {
    if (flow_law->flow(s, E, p, gs[i]) != ref) {
      return true;
    }
  }
  return false;
}

// Rather than make this part of the base class, we just check at some reference values.
bool FlowLawIsPatersonBuddCold(FlowLaw *flow_law, const Config &config,
                               EnthalpyConverter::Ptr EC) {
  static const struct {double s, E, p, gs;} v[] = {
    {1e3, 223, 1e6, 1e-3}, {450000, 475000, 500000, 525000}, {5e4, 268, 5e6, 3e-3}, {1e5, 273, 8e6, 5e-3}};
  PatersonBuddCold cpb("sia_", config, EC); // This is unmodified cold Paterson-Budd
  for (int i=0; i<4; i++) {
    const double left  = flow_law->flow(v[i].s, v[i].E, v[i].p, v[i].gs),
      right =  cpb.flow(v[i].s, v[i].E, v[i].p, v[i].gs);
    if (fabs((left - right)/left)>1.0e-15) {
      return false;
    }
  }
  return true;
}

FlowLaw::FlowLaw(const std::string &prefix, const Config &config,
                 EnthalpyConverter::Ptr EC)
  : m_EC(EC), m_e(1) {

  if (not m_EC) {
    throw RuntimeError("EC is NULL in FlowLaw::FlowLaw()");
  }

  m_standard_gravity   = config.get_double("standard_gravity");
  m_ideal_gas_constant = config.get_double("ideal_gas_constant");

  m_rho                = config.get_double("ice_density");
  m_beta_CC_grad       = config.get_double("beta_CC") * m_rho * m_standard_gravity;
  m_melting_point_temp = config.get_double("water_melting_point_temperature");
  m_e                  = config.get_double(prefix + "enhancement_factor");
  m_n                  = config.get_double(prefix + "Glen_exponent");
  m_viscosity_power    = (1.0 - m_n) / (2.0 * m_n);
  m_hardness_power     = -1.0 / m_n;

  m_A_cold = config.get_double("Paterson_Budd_A_cold");
  m_A_warm = config.get_double("Paterson_Budd_A_warm");
  m_Q_cold = config.get_double("Paterson_Budd_Q_cold");
  m_Q_warm = config.get_double("Paterson_Budd_Q_warm");
  m_crit_temp = config.get_double("Paterson_Budd_critical_temperature");
  m_schoofLen = config.get_double("Schoof_regularizing_length", "m"); // convert to meters
  m_schoofVel = config.get_double("Schoof_regularizing_velocity", "m/s"); // convert to m/s
  m_schoofReg = PetscSqr(m_schoofVel/m_schoofLen);
}

FlowLaw::~FlowLaw() {
  // empty
}

std::string FlowLaw::name() const {
  return m_name;
}

EnthalpyConverter::Ptr FlowLaw::EC() const {
  return m_EC;
}

//! Return the softness parameter A(T) for a given temperature T.
/*! This is not a natural part of all FlowLaw instances.   */
double FlowLaw::softness_paterson_budd(double T_pa) const {
  if (T_pa < m_crit_temp) {
    return m_A_cold * exp(-m_Q_cold/(m_ideal_gas_constant * T_pa));
  }
  return m_A_warm * exp(-m_Q_warm/(m_ideal_gas_constant * T_pa));
}

//! The flow law itself.
double FlowLaw::flow(double stress, double enthalpy,
                     double pressure, double gs) const {
  return this->flow_impl(stress, enthalpy, pressure, gs);
}

double FlowLaw::flow_impl(double stress, double enthalpy,
                          double pressure, double /* gs */) const {
  return softness(enthalpy, pressure) * pow(stress, m_n-1);
}

double FlowLaw::softness(double E, double p) const {
  return this->softness_impl(E, p);
}

double FlowLaw::hardness(double E, double p) const {
  return this->hardness_impl(E, p);
}

double FlowLaw::hardness_impl(double E, double p) const {
  return pow(softness(E, p), m_hardness_power);
}

void averaged_hardness_vec(const FlowLaw &ice,
                           const IceModelVec2S &thickness,
                           const IceModelVec3  &enthalpy,
                           IceModelVec2S &result) {

  const IceGrid &grid = *thickness.get_grid();

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(result);
  list.add(enthalpy);

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Evaluate column integrals in flow law at every quadrature point's column
      double H = thickness(i,j);
      const double *enthColumn = enthalpy.get_column(i, j);
      result(i,j) = averaged_hardness(ice, H, grid.kBelowHeight(H),
                                      &(grid.z()[0]), enthColumn);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  result.update_ghosts();
}


//! Computes vertical average of `B(E, p)` ice hardness, namely @f$\bar B(E, p)@f$.
/*!
 * See comment for hardness(). Note `E[0], ..., E[kbelowH]` must be valid.
 */
double averaged_hardness(const FlowLaw &ice,
                         double thickness, int kbelowH,
                         const double *zlevels,
                         const double *enthalpy) {
  double B = 0;

  EnthalpyConverter &EC = *ice.EC();

  // Use trapezoidal rule to integrate from 0 to zlevels[kbelowH]:
  if (kbelowH > 0) {
    double
      p0 = EC.pressure(thickness),
      E0 = enthalpy[0],
      h0 = ice.hardness(E0, p0); // ice hardness at the left endpoint

    for (int i = 1; i <= kbelowH; ++i) { // note the "1" and the "<="
      const double
        p1 = EC.pressure(thickness - zlevels[i]), // pressure at the right endpoint
        E1 = enthalpy[i], // enthalpy at the right endpoint
        h1 = ice.hardness(E1, p1); // ice hardness at the right endpoint

      // The trapezoid rule sans the "1/2":
      B += (zlevels[i] - zlevels[i-1]) * (h0 + h1);

      h0 = h1;
    }
  }

  // Add the "1/2":
  B *= 0.5;

  // use the "rectangle method" to integrate from
  // zlevels[kbelowH] to thickness:
  double
    depth = thickness - zlevels[kbelowH],
    p = EC.pressure(depth);

  B += depth * ice.hardness(enthalpy[kbelowH], p);

  // Now B is an integral of ice hardness; next, compute the average:
  if (thickness > 0) {
    B = B / thickness;
  } else {
    B = 0;
  }

  return B;
}

/*!
  This constructor just sets flow law factor for nonzero water content, from
  \ref AschwandenBlatter and \ref LliboutryDuval1985.
*/
GPBLD::GPBLD(const std::string &prefix,
             const Config &config, EnthalpyConverter::Ptr EC)
  : FlowLaw(prefix, config, EC) {
  m_name = "Glen-Paterson-Budd-Lliboutry-Duval";

  m_T_0              = config.get_double("water_melting_point_temperature");    // K
  m_water_frac_coeff = config.get_double("gpbld_water_frac_coeff");
  m_water_frac_observed_limit
    = config.get_double("gpbld_water_frac_observed_limit");
}

//! The softness factor in the Glen-Paterson-Budd-Lliboutry-Duval flow law.  For constitutive law form.
/*!
  This is a modification of Glen-Paterson-Budd ice, which is PatersonBudd.  In particular, if
  \f$A()\f$ is the softness factor for PatersonBudd, if \f$E\f$ is the enthalpy, and \f$p\f$ is
  the pressure then the softness we compute is
  \f[A = A(T_{pa}(E, p))(1+184\omega).\f]
  The pressure-melting temperature \f$T_{pa}(E, p)\f$ is computed by pressure_adjusted_temperature().
*/
double GPBLD::softness_impl(double enthalpy, double pressure) const {
  const double E_s = m_EC->enthalpy_cts(pressure);
  if (enthalpy < E_s) {       // cold ice
    double T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure);
    return softness_paterson_budd(T_pa);
  } else { // temperate ice
    double omega = m_EC->water_fraction(enthalpy, pressure);
    // as stated in \ref AschwandenBuelerBlatter, cap omega at max of observations:
    omega = std::min(omega, m_water_frac_observed_limit);
    // next line implements eqn (23) in \ref AschwandenBlatter2009
    return softness_paterson_budd(m_T_0) * (1.0 + m_water_frac_coeff * omega);
  }
}

// PatersonBudd

PatersonBudd::PatersonBudd(const std::string &prefix,
                           const Config &config,
                           EnthalpyConverter::Ptr EC)
  : FlowLaw(prefix, config, EC) {
  m_name = "Paterson-Budd";
}

PatersonBudd::~PatersonBudd() {
  // empty
}

/*! Converts enthalpy to temperature and uses the Paterson-Budd formula. */
double PatersonBudd::softness_impl(double E, double pressure) const {
  double T_pa = m_EC->pressure_adjusted_temperature(E, pressure);
  return softness_from_temp(T_pa);
}

/*! Converts enthalpy to temperature and calls flow_from_temp. */
double PatersonBudd::flow_impl(double stress, double E,
                               double pressure, double gs) const {
  double temp = m_EC->temperature(E, pressure);
  return flow_from_temp(stress, temp, pressure, gs);
}

//! The flow law (temperature-dependent version).
double PatersonBudd::flow_from_temp(double stress, double temp,
                                    double pressure, double /*gs*/) const {
  // pressure-adjusted temperature:
  const double T_pa = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  return softness_from_temp(T_pa) * pow(stress, m_n-1);
}

PatersonBuddCold::PatersonBuddCold(const std::string &prefix,
                   const Config &config,
                   EnthalpyConverter::Ptr EC)
  : PatersonBudd(prefix, config, EC) {
  m_name = "Paterson-Budd (cold case)";
}

double PatersonBuddCold::tempFromSoftness(double myA) const {
  return - Q() / (m_ideal_gas_constant * (log(myA) - log(A())));
}

PatersonBuddCold::~PatersonBuddCold() {
  // empty
}

PatersonBuddWarm::PatersonBuddWarm(const std::string &prefix,
                   const Config &config, EnthalpyConverter::Ptr EC)
  : PatersonBuddCold(prefix, config, EC) {
  m_name = "Paterson-Budd (warm case)";
}

PatersonBuddWarm::~PatersonBuddWarm() {
}


// IsothermalGlen

IsothermalGlen::IsothermalGlen(const std::string &prefix,
                               const Config &config, EnthalpyConverter::Ptr EC)
  : PatersonBudd(prefix, config, EC) {
  m_name = "isothermal Glen";
  
  m_softness_A = config.get_double("ice_softness");
  m_hardness_B = pow(m_softness_A, m_hardness_power);
}

// Hooke

Hooke::Hooke(const std::string &prefix,
             const Config &config, EnthalpyConverter::Ptr EC)
  : PatersonBudd(prefix, config, EC) {
  m_name = "Hooke";

  m_Q_Hooke  = config.get_double("Hooke_Q");
  m_A_Hooke  = config.get_double("Hooke_A");
  m_C_Hooke  = config.get_double("Hooke_C");
  m_K_Hooke  = config.get_double("Hooke_k");
  m_Tr_Hooke = config.get_double("Hooke_Tr");
}

Hooke::~Hooke() {
  // empty
}

double Hooke::softness_from_temp(double T_pa) const {
  return m_A_Hooke * exp(-m_Q_Hooke/(m_ideal_gas_constant * T_pa)
                         + 3.0 * m_C_Hooke * pow(m_Tr_Hooke - T_pa, -m_K_Hooke));
}

// Goldsby-Kohlstedt (forward) ice flow law

GoldsbyKohlstedt::GoldsbyKohlstedt(const std::string &prefix,
                                   const Config &config, EnthalpyConverter::Ptr EC)
  : FlowLaw(prefix, config, EC) {
  m_name = "Goldsby-Kohlstedt / Paterson-Budd (hybrid)";

  m_V_act_vol      = -13.e-6;   // m^3/mol
  m_d_grain_size   = 1.0e-3;    // m  (see p. ?? of G&K paper)
  //--- dislocation creep ---
  m_disl_crit_temp = 258.0;     // Kelvin
  //disl_A_cold    = 4.0e5;                  // MPa^{-4.0} s^{-1}
  //disl_A_warm    = 6.0e28;                 // MPa^{-4.0} s^{-1}
  m_disl_A_cold    = 4.0e-19;   // Pa^{-4.0} s^{-1}
  m_disl_A_warm    = 6.0e4;     // Pa^{-4.0} s^{-1} (GK)
  m_disl_n         = 4.0;       // stress exponent
  m_disl_Q_cold    = 60.e3;     // J/mol Activation energy
  m_disl_Q_warm    = 180.e3;    // J/mol Activation energy (GK)
  //--- grain boundary sliding ---
  m_gbs_crit_temp  = 255.0;     // Kelvin
  //gbs_A_cold     = 3.9e-3;                  // MPa^{-1.8} m^{1.4} s^{-1}
  //gbs_A_warm     = 3.e26;                   // MPa^{-1.8} m^{1.4} s^{-1}
  m_gbs_A_cold     = 6.1811e-14; // Pa^{-1.8} m^{1.4} s^{-1}
  m_gbs_A_warm     = 4.7547e15; // Pa^{-1.8} m^{1.4} s^{-1}
  m_gbs_n          = 1.8;       // stress exponent
  m_gbs_Q_cold     = 49.e3;     // J/mol Activation energy
  m_gbs_Q_warm     = 192.e3;    // J/mol Activation energy
  m_p_grain_sz_exp = 1.4;       // from Peltier
  //--- easy slip (basal) ---
  //basal_A        = 5.5e7;                      // MPa^{-2.4} s^{-1}
  m_basal_A        = 2.1896e-7; // Pa^{-2.4} s^{-1}
  m_basal_n        = 2.4;       // stress exponent
  m_basal_Q        = 60.e3;     // J/mol Activation energy
  //--- diffusional flow ---
  m_diff_crit_temp = 258.0;     // when to use enhancement factor
  m_diff_V_m       = 1.97e-5;   // Molar volume (m^3/mol)
  m_diff_D_0v      = 9.10e-4;   // Preexponential volume diffusion (m^2/s)
  m_diff_Q_v       = 59.4e3;    // activation energy, vol. diff. (J/mol)
  m_diff_D_0b      = 5.8e-4;    // preexponential grain boundary coeff.
  m_diff_Q_b       = 49.e3;     // activation energy, g.b. (J/mol)
  m_diff_delta     = 9.04e-10;  // grain boundary width (m)
}

double GoldsbyKohlstedt::flow_impl(double stress, double E,
                                   double pressure, double grainsize) const {
  double temp = m_EC->temperature(E, pressure);
  return flow_from_temp(stress, temp, pressure, grainsize);
}

double GoldsbyKohlstedt::averaged_hardness_impl(double, int,
                                                const double *,
                                                const double *) const {

  throw std::runtime_error("double GoldsbyKohlstedt::averaged_hardness is not implemented");

#ifndef __GNUC__
  return 0;
#endif
}

double GoldsbyKohlstedt::hardness_impl(double enthalpy, double pressure) const {

  // We use the Paterson-Budd relation for the hardness parameter. It would be nice if we didn't
  // have to, but we currently need ice hardness to compute the strain heating. See
  // SIAFD::compute_volumetric_strain_heating().
  double
    T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure),
    softness = softness_paterson_budd(T_pa);

  return pow(softness, m_hardness_power);
}

double GoldsbyKohlstedt::softness_impl(double , double) const {
  throw std::runtime_error("double GoldsbyKohlstedt::softness is not implemented");

#ifndef __GNUC__
  return 0;
#endif
}

/*!
  This is the (forward) Goldsby-Kohlstedt flow law.  See:
  D. L. Goldsby & D. L. Kohlstedt (2001), "Superplastic deformation
  of ice: experimental observations", J. Geophys. Res. 106(M6), 11017-11030.
*/
double GoldsbyKohlstedt::flow_from_temp(double stress, double temp,
                                        double pressure, double gs) const {
  double eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;

  if (fabs(stress) < 1e-10) {
    return 0;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double pV = pressure * m_V_act_vol;
  const double RT = m_ideal_gas_constant * T;
  // Diffusional Flow
  const double diff_D_v = m_diff_D_0v * exp(-m_diff_Q_v/RT);
  diff_D_b = m_diff_D_0b * exp(-m_diff_Q_b/RT);
  if (T > m_diff_crit_temp) {
    diff_D_b *= 1000; // Coble creep scaling
  }
  eps_diff = 14 * m_diff_V_m *
    (diff_D_v + M_PI * m_diff_delta * diff_D_b / gs) / (RT*PetscSqr(gs));
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_warm + pV)/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_cold + pV)/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-(m_basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_warm + pV)/RT);
  } else {
    eps_gbs = m_gbs_A_cold * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_cold + pV)/RT);
  }

  return eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}


/*****************
THE NEXT PROCEDURE REPEATS CODE; INTENDED ONLY FOR DEBUGGING
*****************/
GKparts GoldsbyKohlstedt::flowParts(double stress, double temp, double pressure) const {
  double gs, eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;
  GKparts p;

  if (fabs(stress) < 1e-10) {
    p.eps_total=0.0;
    p.eps_diff=0.0; p.eps_disl=0.0; p.eps_gbs=0.0; p.eps_basal=0.0;
    return p;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double pV = pressure * m_V_act_vol;
  const double RT = m_ideal_gas_constant * T;
  // Diffusional Flow
  const double diff_D_v = m_diff_D_0v * exp(-m_diff_Q_v/RT);
  diff_D_b = m_diff_D_0b * exp(-m_diff_Q_b/RT);
  if (T > m_diff_crit_temp) {
    diff_D_b *= 1000; // Coble creep scaling
  }
  gs = m_d_grain_size;
  eps_diff = 14 * m_diff_V_m *
    (diff_D_v + M_PI * m_diff_delta * diff_D_b / gs) / (RT*PetscSqr(gs));
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_warm + pV)/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_cold + pV)/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-(m_basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_warm + pV)/RT);
  } else {
    eps_gbs = m_gbs_A_cold * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_cold + pV)/RT);
  }

  p.eps_diff=eps_diff;
  p.eps_disl=eps_disl;
  p.eps_basal=eps_basal;
  p.eps_gbs=eps_gbs;
  p.eps_total=eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
  return p;
}
/*****************/

GoldsbyKohlstedtStripped::GoldsbyKohlstedtStripped(const std::string &prefix,
                                                   const Config &config,
                                                   EnthalpyConverter::Ptr EC)
  : GoldsbyKohlstedt(prefix, config, EC) {
  m_name = "Goldsby-Kohlstedt / Paterson-Budd (hybrid, simplified)";

  m_d_grain_size_stripped = 3.0e-3;  // m; = 3mm  (see Peltier et al 2000 paper)
}


double GoldsbyKohlstedtStripped::flow_from_temp(double stress, double temp, double pressure, double) const {
  // note value of gs is ignored
  // note pressure only effects the temperature; the "P V" term is dropped
  // note no diffusional flow
  double eps_disl, eps_basal, eps_gbs;

  if (fabs(stress) < 1e-10) {
    return 0;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double RT = m_ideal_gas_constant * T;
  // NO Diffusional Flow
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-m_disl_Q_warm/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-m_disl_Q_cold/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-m_basal_Q/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm *
      (pow(stress, m_gbs_n-1) / pow(m_d_grain_size_stripped, m_p_grain_sz_exp)) *
      exp(-m_gbs_Q_warm/RT);
  } else {
    eps_gbs = m_gbs_A_cold *
      (pow(stress, m_gbs_n-1) / pow(m_d_grain_size_stripped, m_p_grain_sz_exp)) *
      exp(-m_gbs_Q_cold/RT);
  }

  return eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}

//! \brief Computes the regularized effective viscosity and its derivative with respect to the
//! second invariant \f$ \gamma \f$.
/*!
 *
 * @f{align*}{
 * \nu &= \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)},\\
 * \diff{\nu}{\gamma} &= \frac{1}{2} B \cdot \frac{1-n}{2n} \cdot \left(\epsilon + \gamma \right)^{(1-n)/(2n) - 1}, \\
 * &= \frac{1-n}{2n} \cdot \frac{1}{2} B \left( \epsilon + \gamma \right)^{(1-n)/(2n)} \cdot \frac{1}{\epsilon + \gamma}, \\
 * &= \frac{1-n}{2n} \cdot \frac{\nu}{\epsilon + \gamma}.
 * @f}
 * Here @f$ \gamma @f$ is the second invariant
 * @f{align*}{
 * \gamma &= \frac{1}{2} D_{ij} D_{ij}\\
 * &= \frac{1}{2}\, ((u_x)^2 + (v_y)^2 + (u_x + v_y)^2 + \frac{1}{2}\, (u_y + v_x)^2) \\
 * @f}
 * and
 * @f[ D_{ij}(\mathbf{u}) = \frac{1}{2}\left(\diff{u_{i}}{x_{j}} + \diff{u_{j}}{x_{i}}\right). @f]
 *
 * Either one of \c nu and \c dnu can be NULL if the corresponding output is not needed.
 *
 * \param[in] hardness ice hardness
 * \param[in] gamma the second invariant
 * \param[out] nu effective viscosity
 * \param[out] dnu derivative of \f$ \nu \f$ with respect to \f$ \gamma \f$
 */
void FlowLaw::effective_viscosity(double hardness, double gamma,
                                         double *nu, double *dnu) const {
  const double
    my_nu = 0.5 * hardness * pow(m_schoofReg + gamma, m_viscosity_power);

  if (PetscLikely(nu != NULL)) {
    *nu = my_nu;
  }

  if (PetscLikely(dnu != NULL)) {
    *dnu = m_viscosity_power * my_nu / (m_schoofReg + gamma);
  }
}

double FlowLaw::exponent() const {
  return m_n;
}

double FlowLaw::enhancement_factor() const {
  return m_e;
}

} // end of namespace rheology
} // end of namespace pism
