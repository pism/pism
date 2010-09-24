// Copyright (C) 2004-2010 Jed Brown, Ed Bueler, and Constantine Khroulev
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
  IceFlowLaw(MPI_Comm c,const char pre[], const NCConfigVariable &config);
  virtual ~IceFlowLaw();
  virtual PetscErrorCode setFromOptions();

  virtual PetscReal effectiveViscosity(PetscReal hardness,
                                       PetscReal u_x, PetscReal u_y,
                                       PetscReal v_x, PetscReal v_y) const;
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
    k,            //!< thermal conductivity
    c_p,          //!< specific heat capacity
    latentHeat,   //!< latent heat capacity
    meltingTemp;  //!< melting temperature
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
  EnthalpyConverter *EC;
};

// Helper functions:
PetscTruth IceFlowLawIsPatersonBuddCold(IceFlowLaw *, const NCConfigVariable &);
PetscTruth IceFlowLawUsesGrainSize(IceFlowLaw *);

//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
  See [\ref AschwandenBlatter]. The basic references are [\ref Glen] and [\ref
  PatersonBudd] and [\ref LliboutryDuval1985].
 */
class GPBLDIce : public IceFlowLaw {
public:
  GPBLDIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual ~GPBLDIce() {}

  virtual PetscErrorCode setFromOptions();
  virtual PetscReal softnessParameter_from_enth(PetscReal enthalpy,
                                                PetscReal pressure) const;
protected:
  PetscReal T_0, water_frac_coeff;
};

//! Derived class of IceFlowLaw for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceFlowLaw {
public:
  ThermoGlenIce(MPI_Comm c, const char pre[], const NCConfigVariable &config)
    : IceFlowLaw(c, pre, config) {}
  virtual ~ThermoGlenIce() {}

  // This also takes care of hardnessParameter_from_enth
  virtual PetscReal softnessParameter_from_enth(PetscReal enthalpy, PetscReal pressure) const;

  // expose Paterson-Budd softness and hardness formulas
  virtual PetscReal softnessParameter(PetscReal T_pa) const
  { return softnessParameter_paterson_budd(T_pa); }
  virtual PetscReal hardnessParameter(PetscReal T_pa) const
  { return pow(softnessParameter(T_pa), -1.0/n); }

  virtual PetscReal flow_from_enth(PetscReal stress, PetscReal E,
                                   PetscReal pressure, PetscReal gs) const;

  // create a special temperature-dependent method
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal pressure, PetscReal gs) const;
};

//! Isothermal Glen ice allowing extra customization.
class CustomGlenIce : public ThermoGlenIce {
public:
  CustomGlenIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual ~CustomGlenIce() {}

  virtual PetscReal softnessParameter_from_enth(PetscReal, PetscReal) const
  { return softness_A; }
  virtual PetscReal softnessParameter(PetscReal) const
  { return softness_A; }

  virtual PetscReal hardnessParameter_from_enth(PetscReal, PetscReal) const
  { return hardness_B; }
  virtual PetscReal hardnessParameter(PetscReal) const
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
  virtual PetscErrorCode setFromOptions();

protected:
  PetscReal softness_A, hardness_B;
};

//! The Hooke flow law.
class HookeIce : public ThermoGlenIce {
public:
  HookeIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual ~HookeIce() {}
  virtual PetscReal softnessParameter(PetscReal T_pa) const;
protected:
  PetscReal A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke; // constants from Hooke (1981)
  // R_Hooke is the ideal_gas_constant.
};

class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c, const char pre[], const NCConfigVariable &config)
    : ThermoGlenIce(c, pre, config) {}
  virtual ~ThermoGlenArrIce() {}
  virtual PetscReal A() const { return A_cold; }
  virtual PetscReal Q() const { return Q_cold; }

  // takes care of hardnessParameter...
  virtual PetscReal softnessParameter(PetscReal T_pa) const
  { return A() * exp(-Q()/(ideal_gas_constant * T_pa)); }

  //! Return the temperature T corresponding to a given value A=A(T).
  PetscReal tempFromSoftness(PetscScalar myA) const
  { return - Q() / (ideal_gas_constant * (log(myA) - log(A()))); }

  // ignores pressure and uses non-pressure-adjusted temperature
  virtual PetscReal flow_from_temp(PetscReal stress, PetscReal temp,
                                   PetscReal , PetscReal ) const
  { return softnessParameter(temp) * pow(stress,n-1); }
};

class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
  ThermoGlenArrIceWarm(MPI_Comm c, const char pre[], const NCConfigVariable &config)
    : ThermoGlenArrIce(c, pre, config) {}
  virtual ~ThermoGlenArrIceWarm() {}

  virtual PetscReal A() const { return A_warm; }
  virtual PetscReal Q() const { return Q_warm; }
};

#endif // __flowlaws_hh
