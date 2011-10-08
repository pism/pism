// Copyright (C) 2004-2011 Jed Brown, Ed Bueler, and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include <petsc.h>
#include "enthalpyConverter.hh"
#include "NCVariable.hh"

// This uses the definition of second invariant from Hutter and several others, namely
// \f$ \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline PetscReal secondInvariant(PetscReal u_x, PetscReal u_y,
                                        PetscReal v_x, PetscReal v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline PetscReal secondInvariantDu(const PetscReal Du[])
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

  \note IceFlowLaw derived classes should implement hardnessParameter... in
  terms of softnessParameter... That way in many cases we only need to
  re-implement softnessParameter... to turn one flow law into another.
 */
class IceFlowLaw {
public:
  IceFlowLaw(MPI_Comm c,const char pre[], const NCConfigVariable &config,
             EnthalpyConverter *EC);
  virtual ~IceFlowLaw() {}
  virtual PetscErrorCode setFromOptions();

  virtual PetscReal effectiveViscosity(PetscReal hardness,
                                       PetscReal u_x, PetscReal u_y,
                                       PetscReal v_x, PetscReal v_y) const;

  virtual void effectiveViscosity_with_derivative(PetscReal hardness, const PetscReal Du[],
                                                  PetscReal *nu, PetscReal *dnu) const;

  virtual PetscReal averagedHardness_from_enth(PetscReal thickness,
                                               PetscInt kbelowH,
                                               const PetscReal *zlevels,
                                               const PetscReal *enthalpy) const;
  virtual PetscReal exponent() const { return n; }


  virtual PetscReal hardnessParameter_from_enth(PetscReal E, PetscReal p) const;
  virtual PetscReal softnessParameter_from_enth(PetscReal E, PetscReal p) const = 0;
  virtual PetscReal flow_from_enth(PetscReal stress, PetscReal E,
                                   PetscReal pressure, PetscReal grainsize) const;

  // ideally these would be protected or private:
  PetscReal rho,          //!< ice density
    beta_CC_grad, //!< Clausius-Clapeyron gradient
    melting_point_temp;  //!< for water, 273.15 K
  EnthalpyConverter *EC;
protected:
  PetscReal softnessParameter_paterson_budd(PetscReal T_pa) const;

  PetscReal schoofLen,schoofVel,schoofReg,
    A_cold, A_warm, Q_cold, Q_warm,  // see Paterson & Budd (1982)
    crit_temp;

  PetscReal standard_gravity,
    ideal_gas_constant,
    n;                          // power law exponent

  MPI_Comm comm;
  char prefix[256];
};

// Helper functions:
PetscTruth IceFlowLawIsPatersonBuddCold(IceFlowLaw *, const NCConfigVariable &, EnthalpyConverter*);
PetscTruth IceFlowLawUsesGrainSize(IceFlowLaw *);

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
 */
class GPBLDIce : public IceFlowLaw {
public:
  GPBLDIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
           EnthalpyConverter *EC);
  virtual ~GPBLDIce() {}

  virtual PetscErrorCode setFromOptions();
  virtual PetscReal softnessParameter_from_enth(PetscReal enthalpy,
                                                PetscReal pressure) const;
protected:
  PetscReal T_0, water_frac_coeff, water_frac_observed_limit;
};

//! Derived class of IceFlowLaw for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceFlowLaw {
public:
  ThermoGlenIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
                EnthalpyConverter *my_EC)
    : IceFlowLaw(c, pre, config, my_EC) {
    n = 3;    // Paterson-Budd has the fixed Glen exponent, so it's hard-wired.
  }
  virtual ~ThermoGlenIce() {}

  // This also takes care of hardnessParameter_from_enth
  virtual PetscReal softnessParameter_from_enth(PetscReal enthalpy, PetscReal pressure) const;

  // expose Paterson-Budd softness and hardness formulas
  virtual PetscReal softnessParameter_from_temp(PetscReal T_pa) const
  { return softnessParameter_paterson_budd(T_pa); }
  virtual PetscReal hardnessParameter_from_temp(PetscReal T_pa) const
  { return pow(softnessParameter_from_temp(T_pa), -1.0/n); }

  virtual PetscReal flow_from_enth(PetscReal stress, PetscReal E,
                                   PetscReal pressure, PetscReal gs) const;

  // create a special temperature-dependent method
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal pressure, PetscReal gs) const;
};

//! Isothermal Glen ice allowing extra customization.
class CustomGlenIce : public ThermoGlenIce {
public:
  CustomGlenIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
                EnthalpyConverter *my_EC);
  virtual ~CustomGlenIce() {}

  virtual PetscReal softnessParameter_from_enth(PetscReal, PetscReal) const
  { return softness_A; }
  virtual PetscReal softnessParameter_from_temp(PetscReal) const
  { return softness_A; }

  virtual PetscReal hardnessParameter_from_enth(PetscReal, PetscReal) const
  { return hardness_B; }
  virtual PetscReal hardnessParameter_from_temp(PetscReal) const
  { return hardness_B; }
  virtual PetscReal averagedHardness_from_enth(PetscReal, PetscInt,
                                               const PetscReal*, const PetscReal*) const
  { return hardness_B; }

  virtual PetscReal flow_from_enth(PetscReal stress, PetscReal,
                                   PetscReal, PetscReal ) const
  { return softness_A * pow(stress,n-1); }

  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal,
                                   PetscReal, PetscReal ) const
  { return softness_A * pow(stress,n-1); }

  virtual void setHardness(PetscReal hardness);
  virtual void setSoftness(PetscReal softness);
  virtual void setExponent(PetscReal exponent);
  virtual void setDensity(PetscReal density);
  virtual void setSchoofRegularization(PetscReal vel_peryear, PetscReal len_km);
protected:
  PetscReal softness_A, hardness_B;
};

//! The Hooke flow law.
class HookeIce : public ThermoGlenIce {
public:
  HookeIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
           EnthalpyConverter *EC);
  virtual ~HookeIce() {}
  virtual PetscReal softnessParameter_from_temp(PetscReal T_pa) const;
protected:
  PetscReal A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

//! Cold case of Paterson-Budd
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
                   EnthalpyConverter *my_EC)
    : ThermoGlenIce(c, pre, config, my_EC) {}
  virtual ~ThermoGlenArrIce() {}
  virtual PetscReal A() const { return A_cold; }
  virtual PetscReal Q() const { return Q_cold; }

  // takes care of hardnessParameter...
  virtual PetscReal softnessParameter_from_temp(PetscReal T_pa) const
  { return A() * exp(-Q()/(ideal_gas_constant * T_pa)); }

  //! Return the temperature T corresponding to a given value A=A(T).
  PetscReal tempFromSoftness(PetscReal myA) const
  { return - Q() / (ideal_gas_constant * (log(myA) - log(A()))); }

  // ignores pressure and uses non-pressure-adjusted temperature
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal , PetscReal ) const
  { return softnessParameter_from_temp(temp) * pow(stress,n-1); }
};

//! Warm case of Paterson-Budd
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  ThermoGlenArrIceWarm(MPI_Comm c, const char pre[],
                       const NCConfigVariable &config, EnthalpyConverter *my_EC)
    : ThermoGlenArrIce(c, pre, config, my_EC) {}
  virtual ~ThermoGlenArrIceWarm() {}

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
inverted-and-vertically-integrated flow law ("effectiveViscosity()").  Only the
former form of the flow law is known for Goldsby-Kohlstedt.  If one can
invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
derived class.
 */
class HybridIce : public ThermoGlenIce {
public:
  HybridIce(MPI_Comm c, const char pre[], const NCConfigVariable &config,
            EnthalpyConverter *my_EC);
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                     PetscReal pressure, PetscReal gs) const;
  GKparts flowParts(PetscReal stress, PetscReal temp, PetscReal pressure) const;

protected:
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

//! Derived class of HybridIce for testing purposes only.
/*! 
HybridIceStripped is a simplification of Goldsby-Kohlstedt.  Compare to that used in 
Peltier et al 2000, which is even simpler.
 */
class HybridIceStripped : public HybridIce {
public:
  HybridIceStripped(MPI_Comm c, const char pre[],
                    const NCConfigVariable &config, EnthalpyConverter *my_EC);
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                     PetscReal pressure, PetscReal gs) const;
protected:
  PetscReal d_grain_size_stripped;
};

#endif // __flowlaws_hh
