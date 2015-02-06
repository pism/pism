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

#ifndef __flowlaws_hh
#define __flowlaws_hh

#include <petscsys.h>
#include "PISMUnits.hh"

namespace pism {

class IceModelVec2S;
class IceModelVec3;

class EnthalpyConverter;
class Config;

// This uses the definition of squared second invariant from Hutter and several others, namely the output is
// \f$ D^2 = \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline double secondInvariant_2D(double u_x, double u_y,
                                        double v_x, double v_y) {
  return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x));
}

// The squared second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline double secondInvariantDu_2D(const double Du[]) {
  return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2]));
}

//! Ice flow laws.
namespace rheology {

//! Abstract class containing the constitutive relation for the flow of ice (of
//! the Paterson-Budd type).
/*!
  This is the interface which most of PISM uses for rheology.

  The current implementation of stress-balance computations in PISM restrict
  possible choices of rheologies to ones that

  \li are power laws

  \li allow factoring out a temperature- (or enthalpy-) dependent ice hardness
  factor

  \li can be represented in the viscosity form

  \note FlowLaw derived classes should implement hardness_parameter... in
  terms of softness_parameter... That way in many cases we only need to
  re-implement softness_parameter... to turn one flow law into another.
*/
class FlowLaw {
public:
  FlowLaw(MPI_Comm c, const std::string &prefix,
             const Config &config,
             const EnthalpyConverter *EC);
  virtual ~FlowLaw() {}

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
  inline void effective_viscosity(double hardness, double gamma,
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

  virtual double averaged_hardness(double thickness,
                                   int kbelowH,
                                   const double *zlevels,
                                   const double *enthalpy) const;

  virtual void averaged_hardness_vec(const IceModelVec2S &thickness,
                                     const IceModelVec3& enthalpy,
                                     IceModelVec2S &hardav) const;

  virtual std::string name() const = 0;
  virtual double exponent() const {
    return m_n;
  }
  virtual double enhancement_factor() const {
    return m_e;
  }

  virtual double hardness_parameter(double E, double p) const;
  virtual double softness_parameter(double E, double p) const = 0;
  virtual double flow(double stress, double E,
                      double pressure, double grainsize) const;

protected:
  double m_rho,          //!< ice density
    m_beta_CC_grad, //!< Clausius-Clapeyron gradient
    m_melting_point_temp;  //!< for water, 273.15 K
  const EnthalpyConverter *m_EC;

  double softness_parameter_paterson_budd(double T_pa) const;

  double m_schoofLen, m_schoofVel, m_schoofReg, m_viscosity_power,
    m_hardness_power,
    m_A_cold, m_A_warm, m_Q_cold, m_Q_warm,  // see Paterson & Budd (1982)
    m_crit_temp;

  double m_standard_gravity,
    m_ideal_gas_constant,
    m_e,                          // flow enhancement factor
    m_n;                          // power law exponent

  MPI_Comm m_com;
  std::string m_prefix;           // option (parameter) prefix
};

// Helper functions:
bool FlowLawIsPatersonBuddCold(FlowLaw *, const Config &,
                                  const EnthalpyConverter*);
bool FlowLawUsesGrainSize(FlowLaw *);

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
*/
class GPBLD : public FlowLaw {
public:
  GPBLD(MPI_Comm c, const std::string &prefix,
           const Config &config,
           const EnthalpyConverter *EC);
  virtual ~GPBLD() {}

  virtual double softness_parameter(double enthalpy,
                                    double pressure) const;
  virtual std::string name() const {
    return "Glen-Paterson-Budd-Lliboutry-Duval";
  }
protected:
  double T_0, water_frac_coeff, water_frac_observed_limit;
};

//! Derived class of FlowLaw for Paterson-Budd (1982)-Glen ice.
class PatersonBudd : public FlowLaw {
public:
  PatersonBudd(MPI_Comm c, const std::string &prefix,
                const Config &config,
                const EnthalpyConverter *my_EC)
    : FlowLaw(c, prefix, config, my_EC) {
  }
  virtual ~PatersonBudd() {}

  // This also takes care of hardness_parameter
  virtual double softness_parameter(double enthalpy, double pressure) const;

  virtual double flow(double stress, double E,
                      double pressure, double gs) const;
  virtual std::string name() const {
    return "Paterson-Budd";
  }

protected:
  virtual double softness_parameter_from_temp(double T_pa) const {
    return softness_parameter_paterson_budd(T_pa);
  }

  virtual double hardness_parameter_from_temp(double T_pa) const {
    return pow(softness_parameter_from_temp(T_pa), m_hardness_power);
  }

  // special temperature-dependent method
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;
};

//! Isothermal Glen ice allowing extra customization.
class IsothermalGlen : public PatersonBudd {
public:
  IsothermalGlen(MPI_Comm c, const std::string &prefix,
                    const Config &config,
                    const EnthalpyConverter *my_EC);
  virtual ~IsothermalGlen() {}

  virtual double averaged_hardness(double, int,
                                   const double*, const double*) const {
    return m_hardness_B;
  }

  virtual double flow(double stress, double, double, double) const {
    return m_softness_A * pow(stress, m_n-1);
  }

  virtual double softness_parameter(double, double) const {
    return m_softness_A;
  }

  virtual double hardness_parameter(double, double) const {
    return m_hardness_B;
  }

  virtual std::string name() const {
    return "isothermal Glen";
  }

protected:
  virtual double flow_from_temp(double stress, double, double, double) const {
    return m_softness_A * pow(stress,m_n-1);
  }

protected:
  double m_softness_A, m_hardness_B;
};

//! The Hooke flow law.
class Hooke : public PatersonBudd {
public:
  Hooke(MPI_Comm c, const std::string &prefix,
           const Config &config,
           const EnthalpyConverter *EC);
  virtual ~Hooke() {}
  virtual std::string name() const {
    return "Hooke";
  }
protected:
  virtual double softness_parameter_from_temp(double T_pa) const;

  double m_A_Hooke, m_Q_Hooke, m_C_Hooke, m_K_Hooke, m_Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

//! Cold case of Paterson-Budd
class PatersonBuddCold : public PatersonBudd {
public:
  PatersonBuddCold(MPI_Comm c, const std::string &prefix,
                   const Config &config,
                   const EnthalpyConverter *my_EC)
    : PatersonBudd(c, prefix, config, my_EC) {}
  virtual ~PatersonBuddCold() {}

  //! Return the temperature T corresponding to a given value A=A(T).
  double tempFromSoftness(double myA) const {
    return - Q() / (m_ideal_gas_constant * (log(myA) - log(A())));
  }

  virtual std::string name() const {
    return "Paterson-Budd (cold case)";
  }

protected:
  virtual double A() const {
    return m_A_cold;
  }
  virtual double Q() const {
    return m_Q_cold;
  }

  // takes care of hardness_parameter...
  virtual double softness_parameter_from_temp(double T_pa) const {
    return A() * exp(-Q()/(m_ideal_gas_constant * T_pa));
  }


  // ignores pressure and uses non-pressure-adjusted temperature
  virtual double flow_from_temp(double stress, double temp,
                                double , double) const {
    return softness_parameter_from_temp(temp) * pow(stress,m_n-1);
  }
};

//! Warm case of Paterson-Budd
class PatersonBuddWarm : public PatersonBuddCold {
public:
  PatersonBuddWarm(MPI_Comm c, const std::string &prefix,
                       const Config &config, const EnthalpyConverter *my_EC)
    : PatersonBuddCold(c, prefix, config, my_EC) {}
  virtual ~PatersonBuddWarm() {}

  virtual std::string name() const {
    return "Paterson-Budd (warm case)";
  }

protected:
  virtual double A() const {
    return m_A_warm;
  }
  virtual double Q() const {
    return m_Q_warm;
  }
};

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
  GoldsbyKohlstedt(MPI_Comm c, const std::string &prefix,
                      const Config &config,
                      const EnthalpyConverter *my_EC);

  virtual double flow(double stress, double E,
                      double pressure, double grainsize) const;

  virtual void effective_viscosity(double hardness, double gamma,
                                   double *nu, double *dnu) const __attribute__((noreturn));

  virtual double averaged_hardness(double thickness,
                                   int kbelowH,
                                   const double *zlevels,
                                   const double *enthalpy) const __attribute__((noreturn));

  virtual double hardness_parameter(double E, double p) const;

  virtual double softness_parameter(double E, double p) const __attribute__((noreturn));

  virtual std::string name() const {
    return "Goldsby-Kohlstedt / Paterson-Budd (hybrid)";
  }

protected:
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
  GoldsbyKohlstedtStripped(MPI_Comm c, const std::string &prefix,
                              const Config &config, const EnthalpyConverter *my_EC);
  virtual std::string name() const {
    return "Goldsby-Kohlstedt / Paterson-Budd (hybrid, simplified)";
  }

protected:
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;

  double m_d_grain_size_stripped;
};

} // end of namespace rheology
} // end of namespace pism

#endif // __flowlaws_hh
