// Copyright (C) 2004-2013 Jed Brown, Ed Bueler, and Constantine Khroulev
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

class EnthalpyConverter;
class NCConfigVariable;

// This uses the definition of second invariant from Hutter and several others, namely
// \f$ \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline PetscReal secondInvariant_2D(PetscReal u_x, PetscReal u_y,
					   PetscReal v_x, PetscReal v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline PetscReal secondInvariantDu_2D(const PetscReal Du[])
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
  IceFlowLaw(MPI_Comm c, const char pre[], PISMUnitSystem unit_system,
             const NCConfigVariable &config,
             EnthalpyConverter *EC);
 virtual ~IceFlowLaw() {}
  virtual PetscErrorCode setFromOptions();

  virtual void effective_viscosity(PetscReal hardness, PetscReal gamma,
				   PetscReal *nu, PetscReal *dnu) const;

  virtual PetscReal averaged_hardness(PetscReal thickness,
                                      PetscInt kbelowH,
                                      const PetscReal *zlevels,
                                      const PetscReal *enthalpy) const;
  virtual PetscReal exponent() const { return n; }
  virtual PetscReal enhancement_factor() const { return e; }

  virtual PetscReal hardness_parameter(PetscReal E, PetscReal p) const;
  virtual PetscReal softness_parameter(PetscReal E, PetscReal p) const = 0;
  virtual PetscReal flow(PetscReal stress, PetscReal E,
                         PetscReal pressure, PetscReal grainsize) const;

protected:
  PetscReal rho,          //!< ice density
    beta_CC_grad, //!< Clausius-Clapeyron gradient
    melting_point_temp;  //!< for water, 273.15 K
  EnthalpyConverter *EC;

  PetscReal softness_parameter_paterson_budd(PetscReal T_pa) const;

  PetscReal schoofLen,schoofVel,schoofReg,
    A_cold, A_warm, Q_cold, Q_warm,  // see Paterson & Budd (1982)
    crit_temp;

  PetscReal standard_gravity,
    ideal_gas_constant,
    e,                          // flow enhancement factor
    n;                          // power law exponent

  MPI_Comm com;
  char prefix[256];
};

// Helper functions:
PetscBool IceFlowLawIsPatersonBuddCold(IceFlowLaw *, const NCConfigVariable &,
                                       PISMUnitSystem unit_system, EnthalpyConverter*);
PetscBool IceFlowLawUsesGrainSize(IceFlowLaw *);

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
 */
class GPBLDIce : public IceFlowLaw {
public:
  GPBLDIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
           const NCConfigVariable &config,
           EnthalpyConverter *EC);
  virtual ~GPBLDIce() {}

  virtual PetscErrorCode setFromOptions();
  virtual PetscReal softness_parameter(PetscReal enthalpy,
                                       PetscReal pressure) const;
protected:
  PetscReal T_0, water_frac_coeff, water_frac_observed_limit;
};

//! Derived class of IceFlowLaw for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceFlowLaw {
public:
  ThermoGlenIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
                const NCConfigVariable &config,
                EnthalpyConverter *my_EC)
    : IceFlowLaw(c, pre, s, config, my_EC) {
  }
  virtual ~ThermoGlenIce() {}

  // This also takes care of hardness_parameter
  virtual PetscReal softness_parameter(PetscReal enthalpy, PetscReal pressure) const;

  virtual PetscReal flow(PetscReal stress, PetscReal E,
                         PetscReal pressure, PetscReal gs) const;

protected:
  virtual PetscReal softness_parameter_from_temp(PetscReal T_pa) const
  { return softness_parameter_paterson_budd(T_pa); }

  virtual PetscReal hardness_parameter_from_temp(PetscReal T_pa) const
  { return pow(softness_parameter_from_temp(T_pa), -1.0/n); }

  // special temperature-dependent method
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal pressure, PetscReal gs) const;
};

//! Isothermal Glen ice allowing extra customization.
class IsothermalGlenIce : public ThermoGlenIce {
public:
  IsothermalGlenIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
                    const NCConfigVariable &config,
                    EnthalpyConverter *my_EC);
  virtual ~IsothermalGlenIce() {}

  virtual PetscReal averaged_hardness(PetscReal, PetscInt,
                                      const PetscReal*, const PetscReal*) const
  { return hardness_B; }

  virtual PetscReal flow(PetscReal stress, PetscReal,
                         PetscReal, PetscReal ) const
  { return softness_A * pow(stress, n-1); }

  virtual PetscReal softness_parameter(PetscReal, PetscReal) const
  { return softness_A; }

  virtual PetscReal hardness_parameter(PetscReal, PetscReal) const
  { return hardness_B; }

protected:
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal,
                                   PetscReal, PetscReal ) const
  { return softness_A * pow(stress,n-1); }

protected:
  PetscReal softness_A, hardness_B;
};

//! The Hooke flow law.
class HookeIce : public ThermoGlenIce {
public:
  HookeIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
           const NCConfigVariable &config,
           EnthalpyConverter *EC);
  virtual ~HookeIce() {}
protected:
  virtual PetscReal softness_parameter_from_temp(PetscReal T_pa) const;

  PetscReal A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

//! Cold case of Paterson-Budd
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
                   const NCConfigVariable &config,
                   EnthalpyConverter *my_EC)
    : ThermoGlenIce(c, pre, s, config, my_EC) {}
  virtual ~ThermoGlenArrIce() {}

  //! Return the temperature T corresponding to a given value A=A(T).
  PetscReal tempFromSoftness(PetscReal myA) const
  { return - Q() / (ideal_gas_constant * (log(myA) - log(A()))); }

protected:
  virtual PetscReal A() const { return A_cold; }
  virtual PetscReal Q() const { return Q_cold; }

  // takes care of hardness_parameter...
  virtual PetscReal softness_parameter_from_temp(PetscReal T_pa) const
  { return A() * exp(-Q()/(ideal_gas_constant * T_pa)); }


  // ignores pressure and uses non-pressure-adjusted temperature
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal , PetscReal ) const
  { return softness_parameter_from_temp(temp) * pow(stress,n-1); }
};

//! Warm case of Paterson-Budd
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  ThermoGlenArrIceWarm(MPI_Comm c, const char pre[], PISMUnitSystem s,
                       const NCConfigVariable &config, EnthalpyConverter *my_EC)
    : ThermoGlenArrIce(c, pre, s, config, my_EC) {}
  virtual ~ThermoGlenArrIceWarm() {}

protected:
  virtual PetscReal A() const { return A_warm; }
  virtual PetscReal Q() const { return Q_warm; }
};

// Hybrid (Goldsby-Kohlstedt/Glen) ice flow law

struct GKparts {
  PetscReal eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
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
  GoldsbyKohlstedtIce(MPI_Comm c, const char pre[], PISMUnitSystem s,
                      const NCConfigVariable &config,
                      EnthalpyConverter *my_EC);

  virtual PetscReal flow(PetscReal stress, PetscReal E,
                         PetscReal pressure, PetscReal grainsize) const;

  virtual void effective_viscosity(PetscReal hardness, PetscReal gamma,
				   PetscReal *nu, PetscReal *dnu) const;

  virtual PetscReal averaged_hardness(PetscReal thickness,
                                      PetscInt kbelowH,
                                      const PetscReal *zlevels,
                                      const PetscReal *enthalpy) const;

  virtual PetscReal hardness_parameter(PetscReal E, PetscReal p) const;

  virtual PetscReal softness_parameter(PetscReal E, PetscReal p) const;

protected:
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal pressure, PetscReal gs) const;
  GKparts flowParts(PetscReal stress, PetscReal temp, PetscReal pressure) const;

  PetscReal  V_act_vol,  d_grain_size,
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
  GoldsbyKohlstedtIceStripped(MPI_Comm c, const char pre[], PISMUnitSystem s,
                              const NCConfigVariable &config, EnthalpyConverter *my_EC);
protected:
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal pressure, PetscReal gs) const;

  PetscReal d_grain_size_stripped;
};

#endif // __flowlaws_hh
