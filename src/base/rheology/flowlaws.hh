// Copyright (C) 2004-2014 Jed Brown, Ed Bueler, and Constantine Khroulev
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
class PISMConfig;

// This uses the definition of squared second invariant from Hutter and several others, namely the output is
// \f$ D^2 = \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline double secondInvariant_2D(double u_x, double u_y,
                                        double v_x, double v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The squared second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline double secondInvariantDu_2D(const double Du[])
{ return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2])); }


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

  \note IceFlowLaw derived classes should implement hardness_parameter... in
  terms of softness_parameter... That way in many cases we only need to
  re-implement softness_parameter... to turn one flow law into another.
*/
class IceFlowLaw {
public:
  IceFlowLaw(MPI_Comm c, const char pre[],
             const PISMConfig &config,
             EnthalpyConverter *EC);
  virtual ~IceFlowLaw() {}
  virtual PetscErrorCode setFromOptions();

  //! \brief Computes effective viscosity and its derivative with respect to the
  //! squared second invariant @f$ \gamma @f$.
  /*!
   *
   * @f{align*}{
   * \nu &= \frac12 B (\gamma + \epsilon)^{\frac{1-n}{2n}}.\\
   * \frac{\partial \nu}{\partial \gamma} &= \frac{1-n}{2n} \cdot \frac12 B (\gamma + \epsilon)^{\frac{1-n}{2n}-1}\\
   * &= \frac{1-n}{2n} \cdot \frac12 B (\gamma + \epsilon)^{\frac{1-n}{2n}}\cdot \frac{1}{\gamma + \epsilon}\\
   * &=\nu \cdot \frac{1-n}{2n} \cdot \frac{1}{\gamma + \epsilon}.
   * @f}
   *
   * Either one of \c nu and \c dnu can be NULL if the corresponding output is not needed.
   *
   * \param[in] hardness ice hardness
   * \param[in] gamma the second invariant @f$ \gamma = \frac{1}{2} D_{ij} D_{ij} @f$ if @f$ D_{ij} @f$ is the strain rate tensor
   * \param[out] nu effective viscosity
   * \param[out] dnu derivative of @f$ \nu @f$ with respect to @f$ \gamma @f$
   */
  inline void effective_viscosity(double hardness, double gamma,
                                  double *nu, double *dnu) const {
    const double
      my_nu = 0.5 * hardness * pow(schoofReg + gamma, viscosity_power);

    if (PetscLikely(nu != NULL))
      *nu = my_nu;
    if (PetscLikely(dnu != NULL))
      *dnu = viscosity_power * my_nu / (schoofReg + gamma);
  }

  virtual double averaged_hardness(double thickness,
                                   int kbelowH,
                                   const double *zlevels,
                                   const double *enthalpy) const;

  virtual PetscErrorCode averaged_hardness_vec(IceModelVec2S &thickness,
                                               IceModelVec3& enthalpy, IceModelVec2S &hardav) const;

  virtual std::string name() const = 0;
  virtual double exponent() const { return n; }
  virtual double enhancement_factor() const { return e; }

  virtual double hardness_parameter(double E, double p) const;
  virtual double softness_parameter(double E, double p) const = 0;
  virtual double flow(double stress, double E,
                      double pressure, double grainsize) const;

protected:
  double rho,          //!< ice density
    beta_CC_grad, //!< Clausius-Clapeyron gradient
    melting_point_temp;  //!< for water, 273.15 K
  EnthalpyConverter *EC;

  double softness_parameter_paterson_budd(double T_pa) const;

  double schoofLen,schoofVel,schoofReg, viscosity_power,
    hardness_power,
    A_cold, A_warm, Q_cold, Q_warm,  // see Paterson & Budd (1982)
    crit_temp;

  double standard_gravity,
    ideal_gas_constant,
    e,                          // flow enhancement factor
    n;                          // power law exponent

  MPI_Comm com;
  char prefix[256];
};

// Helper functions:
PetscBool IceFlowLawIsPatersonBuddCold(IceFlowLaw *, const PISMConfig &,
                                       EnthalpyConverter*);
PetscBool IceFlowLawUsesGrainSize(IceFlowLaw *);

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
*/
class GPBLDIce : public IceFlowLaw {
public:
  GPBLDIce(MPI_Comm c, const char pre[],
           const PISMConfig &config,
           EnthalpyConverter *EC);
  virtual ~GPBLDIce() {}

  virtual double softness_parameter(double enthalpy,
                                    double pressure) const;
  virtual std::string name() const
  { return "Glen-Paterson-Budd-Lliboutry-Duval"; }
protected:
  double T_0, water_frac_coeff, water_frac_observed_limit;
};

//! Derived class of IceFlowLaw for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceFlowLaw {
public:
  ThermoGlenIce(MPI_Comm c, const char pre[],
                const PISMConfig &config,
                EnthalpyConverter *my_EC)
    : IceFlowLaw(c, pre, config, my_EC) {
  }
  virtual ~ThermoGlenIce() {}

  // This also takes care of hardness_parameter
  virtual double softness_parameter(double enthalpy, double pressure) const;

  virtual double flow(double stress, double E,
                      double pressure, double gs) const;
  virtual std::string name() const
  { return "Paterson-Budd"; }

protected:
  virtual double softness_parameter_from_temp(double T_pa) const
  { return softness_parameter_paterson_budd(T_pa); }

  virtual double hardness_parameter_from_temp(double T_pa) const
  { return pow(softness_parameter_from_temp(T_pa), hardness_power); }

  // special temperature-dependent method
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;
};

//! Isothermal Glen ice allowing extra customization.
class IsothermalGlenIce : public ThermoGlenIce {
public:
  IsothermalGlenIce(MPI_Comm c, const char pre[],
                    const PISMConfig &config,
                    EnthalpyConverter *my_EC);
  virtual ~IsothermalGlenIce() {}

  virtual double averaged_hardness(double, int,
                                   const double*, const double*) const
  { return hardness_B; }

  virtual double flow(double stress, double,
                      double, double ) const
  { return softness_A * pow(stress, n-1); }

  virtual double softness_parameter(double, double) const
  { return softness_A; }

  virtual double hardness_parameter(double, double) const
  { return hardness_B; }

  virtual std::string name() const
  { return "isothermal Glen"; }

protected:
  virtual double flow_from_temp(double stress, double,
                                double, double ) const
  { return softness_A * pow(stress,n-1); }

protected:
  double softness_A, hardness_B;
};

//! The Hooke flow law.
class HookeIce : public ThermoGlenIce {
public:
  HookeIce(MPI_Comm c, const char pre[],
           const PISMConfig &config,
           EnthalpyConverter *EC);
  virtual ~HookeIce() {}
  virtual std::string name() const
  { return "Hooke"; }
protected:
  virtual double softness_parameter_from_temp(double T_pa) const;

  double A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

//! Cold case of Paterson-Budd
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c, const char pre[],
                   const PISMConfig &config,
                   EnthalpyConverter *my_EC)
    : ThermoGlenIce(c, pre, config, my_EC) {}
  virtual ~ThermoGlenArrIce() {}

  //! Return the temperature T corresponding to a given value A=A(T).
  double tempFromSoftness(double myA) const
  { return - Q() / (ideal_gas_constant * (log(myA) - log(A()))); }

  virtual std::string name() const
  { return "Paterson-Budd (cold case)"; }

protected:
  virtual double A() const { return A_cold; }
  virtual double Q() const { return Q_cold; }

  // takes care of hardness_parameter...
  virtual double softness_parameter_from_temp(double T_pa) const
  { return A() * exp(-Q()/(ideal_gas_constant * T_pa)); }


  // ignores pressure and uses non-pressure-adjusted temperature
  virtual double flow_from_temp(double stress, double temp,
                                double , double ) const
  { return softness_parameter_from_temp(temp) * pow(stress,n-1); }
};

//! Warm case of Paterson-Budd
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  ThermoGlenArrIceWarm(MPI_Comm c, const char pre[],
                       const PISMConfig &config, EnthalpyConverter *my_EC)
    : ThermoGlenArrIce(c, pre, config, my_EC) {}
  virtual ~ThermoGlenArrIceWarm() {}

  virtual std::string name() const
  { return "Paterson-Budd (warm case)"; }

protected:
  virtual double A() const { return A_warm; }
  virtual double Q() const { return Q_warm; }
};

// Hybrid (Goldsby-Kohlstedt/Glen) ice flow law

struct GKparts {
  double eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};


//! A hybrid of Goldsby-Kohlstedt (2001) ice (constitutive form) and Paterson-Budd (1982)-Glen (viscosity form).
/*!
  Each IceFlowLaw has both a forward flow law in "constitutive law" form ("flow_from_temp()") and an
  inverted-and-vertically-integrated flow law ("effective_viscosity()").  Only the
  former form of the flow law is known for Goldsby-Kohlstedt.  If one can
  invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
  derived class.
*/
class GoldsbyKohlstedtIce : public IceFlowLaw {
public:
  GoldsbyKohlstedtIce(MPI_Comm c, const char pre[],
                      const PISMConfig &config,
                      EnthalpyConverter *my_EC);

  virtual double flow(double stress, double E,
                      double pressure, double grainsize) const;

  virtual void effective_viscosity(double hardness, double gamma,
                                   double *nu, double *dnu) const;

  virtual double averaged_hardness(double thickness,
                                   int kbelowH,
                                   const double *zlevels,
                                   const double *enthalpy) const;

  virtual double hardness_parameter(double E, double p) const;

  virtual double softness_parameter(double E, double p) const;

  virtual std::string name() const
  { return "Goldsby-Kohlstedt / Paterson-Budd (hybrid)"; }

protected:
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;
  GKparts flowParts(double stress, double temp, double pressure) const;

  double  V_act_vol,  d_grain_size,
  //--- diffusional flow ---
    diff_crit_temp, diff_V_m, diff_D_0v, diff_Q_v, diff_D_0b, diff_Q_b, diff_delta,
  //--- dislocation creep ---
    disl_crit_temp, disl_A_cold, disl_A_warm, disl_n, disl_Q_cold, disl_Q_warm,
  //--- easy slip (basal) ---
    basal_A, basal_n, basal_Q,
  //--- grain boundary sliding ---
    gbs_crit_temp, gbs_A_cold, gbs_A_warm, gbs_n, gbs_Q_cold,
    p_grain_sz_exp, gbs_Q_warm;
};

//! Derived class of GoldsbyKohlstedtIce for testing purposes only.
/*!
  GoldsbyKohlstedtIceStripped is a simplification of Goldsby-Kohlstedt. Compare to that
  used in Peltier et al 2000, which is even simpler.
*/
class GoldsbyKohlstedtIceStripped : public GoldsbyKohlstedtIce {
public:
  GoldsbyKohlstedtIceStripped(MPI_Comm c, const char pre[],
                              const PISMConfig &config, EnthalpyConverter *my_EC);
  virtual std::string name() const
  { return "Goldsby-Kohlstedt / Paterson-Budd (hybrid, simplified)"; }

protected:
  virtual double flow_from_temp(double stress, double temp,
                                double pressure, double gs) const;

  double d_grain_size_stripped;
};

} // end of namespace pism

#endif // __flowlaws_hh
