// Copyright (C) 2004-2010 Jed Brown and Ed Bueler
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

#ifndef __materials_hh
#define __materials_hh

#include <petsc.h>
#include "enthalpyConverter.hh"
#include "NCVariable.hh"

/*******************
REGARDING IceFlowLaw and HybridIce:  The main hierarchy is:
  IceFlowLaw <- ThermoGlenIce <- HybridIce <- HybridIceStripped,
where "<-" means "derived class".  IceFlowLaw is a virtual
class; it should never be used "as is".

ThermoGlenIce     means *Paterson-Budd* version of Arhennius relation
ThermoGlenArrIceHooke  means *Hooke* version ...
ThermoGlenArrIce       uses only cold part of Paterson-Budd
ThermoGlenArrIceWarm   uses only warm part of Paterson-Budd
HybridIce         means *Goldsby-Kohlstedt* flow law where vMask=SHEET,
                  and otherwise Paterson-Budd
HybridIceStripped means, where SHEET, G-K without the pressure dependence and without the
                  diffusional part; also grain size fixed at 3mm

Note each IceFlowLaw has both a forward flow law ("flow") and an
inverted-and-vertically-integrated flow law ("effectiveViscosityColumn").  Only the
former form of the flow law is known for Goldsby-Kohlstedt.  If one can
invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
derived class.
*******************/

//! Abstract class containing physical constants and the constitutive relation describing ice.
/*!
This is the interface which most of PISM uses for rheology.
 */
class IceFlowLaw {
public:
  // ideally these would be protected or private:
  PetscScalar rho,          //!< density
              beta_CC_grad, //!< Clausius-Clapeyron gradient
              k,            //!< thermal conductivity
              c_p,          //!< specific heat capacity
              latentHeat,   //!< latent heat capacity
              meltingTemp,  //!< melting temperature
              standard_gravity,
              ideal_gas_constant;

  IceFlowLaw(MPI_Comm c,const char pre[], const NCConfigVariable &config);
  virtual ~IceFlowLaw() {}
  virtual PetscErrorCode setFromOptions() {return 0;}
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const {return 0;}
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const = 0;
  // returns nu * H; it is adapted to a staggered grid so T1,T2 get averaged
  virtual PetscScalar effectiveViscosityColumn(PetscScalar H, PetscInt kbelowH, const PetscScalar *zlevels,
                                               PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y,
                                               const PetscScalar *T1, const PetscScalar *T2) const = 0;
  virtual PetscScalar exponent() const = 0;
  virtual PetscScalar hardnessParameter(PetscScalar T) const = 0;
  virtual PetscScalar averagedHardness(PetscScalar H, PetscInt kbelowH, const PetscScalar zlevels[],
                                      const PetscScalar T[]) const = 0;

protected:
  MPI_Comm comm;
  char prefix[256];
};


PetscTruth IceFlowLawIsPatersonBuddCold(IceFlowLaw *, const NCConfigVariable &);
PetscTruth IceFlowLawUsesGrainSize(IceFlowLaw *);


//! Derived class of IceFlowLaw for which is still isothermal and power law, but has easily setable parameters.
class CustomGlenIce : public IceFlowLaw {
public:
  CustomGlenIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  PetscErrorCode setDensity(PetscReal);
  PetscErrorCode setExponent(PetscReal);
  PetscErrorCode setHardness(PetscReal);
  PetscErrorCode setSoftness(PetscReal);
  PetscErrorCode setSchoofRegularization(PetscReal vel,PetscReal len); // vel has units m/a, len has units km
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar flow(PetscScalar,PetscScalar,PetscScalar,PetscScalar) const;
  virtual PetscScalar effectiveViscosityColumn(PetscScalar,PetscInt,const PetscScalar *,
                                               PetscScalar,PetscScalar,PetscScalar,PetscScalar,
                                               const PetscScalar *,const PetscScalar *) const;
  virtual PetscScalar exponent() const;
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar hardnessParameter(PetscScalar T) const;
  virtual PetscScalar averagedHardness(PetscScalar H, PetscInt kbelowH, const PetscScalar zlevels[],
                                      const PetscScalar T[]) const;
private:
  PetscReal exponent_n,softness_A,hardness_B,schoofVel,schoofLen,schoofReg;
};


//! Derived class of IceFlowLaw for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceFlowLaw {
public:
  ThermoGlenIce(MPI_Comm c,const char pre[], const NCConfigVariable &config);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;
  virtual PetscScalar effectiveViscosityColumn(PetscScalar,PetscInt,const PetscScalar *,
                                               PetscScalar,PetscScalar,PetscScalar,PetscScalar,
                                               const PetscScalar *,const PetscScalar *) const;
  virtual PetscScalar exponent() const;
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar hardnessParameter(PetscScalar T) const;
  virtual PetscScalar averagedHardness(PetscScalar H, PetscInt kbelowH, const PetscScalar zlevels[],
                                      const PetscScalar T[]) const;
protected:
  PetscReal schoofLen,schoofVel,schoofReg,
            A_cold, A_warm, Q_cold, Q_warm,  // these four constants from Paterson & Budd (1982)
            crit_temp, n;
};


//! Glen (1955) and Paterson-Budd (1982) flow law with additional water fraction factor from Lliboutry & Duval (1985).
/*!
See \ref AschwandenBlatter.  The basic references are \ref Glen and \ref PatersonBudd 
and \ref LliboutryDuval1985.
 */
class PolyThermalGPBLDIce : public ThermoGlenIce {
public:
  PolyThermalGPBLDIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual ~PolyThermalGPBLDIce();
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode view(PetscViewer viewer) const;

  /* these are not reimplementations, but new routines.
  to see where these replacements are needed, do in src/base/:
     $ grep ice->flow *.cc
     $ grep ice->effectiveViscosity *.cc
     $ grep ice->softnessParameter *.cc
     $ grep ice->hardnessParameter *.cc
  note for grep: arrow (>) must actually be escaped with backslash
  */
  virtual PetscScalar softnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const;
  virtual PetscScalar hardnessParameterFromEnth(PetscScalar enthalpy, PetscScalar pressure) const;

  virtual PetscScalar flowFromEnth(PetscScalar stress, PetscScalar enthalpy, PetscScalar pressure,
                                       PetscScalar gs) const; // grainsize arg gs not used

  virtual PetscScalar effectiveViscosityColumnFromEnth(
                PetscScalar thickness,  PetscInt kbelowH, const PetscScalar *zlevels,
                PetscScalar u_x,  PetscScalar u_y, PetscScalar v_x,  PetscScalar v_y,
                const PetscScalar *enthalpy1, const PetscScalar *enthalpy2) const;

  // these are used in src/base/ssaJed/ stuff only, so not addressed for now:
  //    integratedStoreSize(), integratedStore(), integratedViscosity()

protected:
  PetscReal T_0, water_frac_coeff;
  EnthalpyConverter *EC;
};


//! Derived class of IceFlowLaw for Hooke (1981)-Glen ice.  Only changes A(T) factor from ThermoGlenIce.
class ThermoGlenIceHooke : public ThermoGlenIce {
public:
  ThermoGlenIceHooke(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual PetscScalar softnessParameter(PetscScalar T) const;
protected:
  PetscReal A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke, R_Hooke; // constants from Hooke (1981)
};


//! Derived class of IceFlowLaw for Arrhenius-Glen ice.  \e Cold case of Paterson-Budd (1982) ice.
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c,const char pre[], const NCConfigVariable &conf) : ThermoGlenIce(c,pre,conf) {}
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar tempFromSoftness(PetscScalar A) const; 
  using ThermoGlenIce::flow;
  virtual PetscScalar flow(PetscScalar,PetscScalar,PetscScalar,PetscScalar) const;
  virtual PetscScalar A() const;  // returns A_cold for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_cold for Paterson-Budd
};


//! Derived class of IceFlowLaw for Arrhenius-Glen ice.  \e Warm case of Paterson-Budd (1982) ice.
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  ThermoGlenArrIceWarm(MPI_Comm c,const char pre[], const NCConfigVariable &conf) : ThermoGlenArrIce(c,pre,conf) {}
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar A() const;  // returns A_warm for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_warm for Paterson-Budd
};



struct GKparts {
  PetscScalar eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};


//! Derived class of IceFlowLaw for a hybrid of Goldsby-Kohlstedt (2001) ice in SIA, with Paterson-Budd (1982)-Glen behavior when needed in viscosity form (e.g. SSA).
class HybridIce : public ThermoGlenIce {
public:
  HybridIce(MPI_Comm c, const char pre[], const NCConfigVariable &config);
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_TRUE; }
  GKparts flowParts(PetscScalar stress, PetscScalar temp, PetscScalar pressure) const;

protected:
  PetscReal  V_act_vol,  d_grain_size,
             //--- diffusional flow ---
             diff_crit_temp, diff_V_m, diff_D_0v, diff_Q_v, diff_D_0b, diff_Q_b, diff_delta,
             //--- dislocation creep ---
             disl_crit_temp, disl_A_cold, disl_A_warm, disl_n, disl_Q_cold, disl_Q_warm,
             //--- easy slip (basal) ---
             basal_A, basal_n, basal_Q,
             //--- grain boundary sliding ---
             gbs_crit_temp, gbs_A_cold, gbs_A_warm, gbs_n, gbs_Q_cold, p_grain_sz_exp, gbs_Q_warm;
};


//! Derived class of HybridIce; for testing only.
/*! 
HybridIceStripped is a simplification of Goldsby-Kohlstedt.  Compare to that used in 
Peltier et al 2000, which is even simpler.
 */
class HybridIceStripped : public HybridIce {
public:
  HybridIceStripped(MPI_Comm c,const char pre[], const NCConfigVariable &config);
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_FALSE; }
protected:
  PetscReal d_grain_size_stripped;
};


//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class IceBasalResistancePlasticLaw {
public:
  IceBasalResistancePlasticLaw(PetscScalar regularizationConstant, bool pseudoPlastic,
                   PetscScalar pseudoExponent, PetscScalar pseudoUThreshold);
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  // Also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy, PetscScalar *drag, PetscScalar *ddrag) const;
  virtual ~IceBasalResistancePlasticLaw() {} // class w virtual methods needs virtual destructor?

  PetscReal   plastic_regularize, pseudo_q, pseudo_u_threshold;
  bool pseudo_plastic;
};


//! Where ice thickness is zero the SSA is no longer "elliptic".  This class provides an extension coefficient to maintain well-posedness/"ellipticity".
/*!
More specifically, the SSA equations are
\latexonly
\def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
\begin{equation*}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  =  - \rho g H \ddx{h},
\end{equation*}
\endlatexonly
and another similar equation.  Schoof \ref SchoofStream shows that these PDEs
are the variational equations for a functional.

The quantity \f$\nu H\f$ is the nonlinear coefficient in this (roughly-speaking)
elliptic pair of PDEs.  Conceptually it is a membrane strength.  Well-posedness of the SSA problem 
requires either a precisely-defined boundary and an appropriate boundary condition
\e or a nonzero value of \f$\nu H\f$ at all points.  This class provides that nonzero value.
 */
class SSAStrengthExtension {
public:
  SSAStrengthExtension();
  virtual ~SSAStrengthExtension();
  //! Set strength with units (viscosity times thickness).
  virtual PetscErrorCode set_notional_strength(PetscReal my_nuH);
  //! Set minimum thickness to trigger use of extension.
  virtual PetscErrorCode set_min_thickness(PetscReal my_min_thickness);
  virtual PetscReal      notional_strength() const;           //!< Returns strength with units (viscosity times thickness).
  virtual PetscReal      min_thickness_for_extension() const; //!< Returns minimum thickness to trigger use of extension.
private:
  PetscReal  min_thickness, nuH;
};



#define ICE_CUSTOM  "custom"        /* Plain isothermal Glen with customizable parameters */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_HYBRID  "hybrid"        /* Goldsby-Kohlstedt for SIA, PB for SSA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

/*!
Excerpt from Jed email to ELB, 2/21 (with fixes to reflect present names of objects (CK)):

IceModel creates an IceFlowLawFactory in it's constructor.  Derived classes
can set options, and driver programs can get a reference and set
options.  If no IceModel::ice object is created before
IceModel::setFromOptions() is called, then IceModel will use the factory
to create ice (whatever options have been set by derived classes serve
as defaults, the user can override these options using -ice_*, run with
-help to see options specific to that particular type).  If the derived
class does not want options to be able to influence ice type (for
instance, to guarantee a particular type), then it can create the ice
before calling IceModel::setFromOptions().  See for instance
IceCompModel and IceExactSSAModel.  This amount of control might not be
desirable in which case just strip it out.

You can check whether an ice type implements a particular interface
(ThermoGlenIce for instance) by using dynamic_cast and checking if the
returned pointer is non-NULL.

You can add ice types without modifying materials.{hh,cc}.  See the
usage in MISMIP, you create a trivial creation function and register it
with the factory.

To test whether an IceFlowLaw object is true to a particular constitutive
relation, it's best to just create a reference object and test at a few
values.  Since I envision IceFlowLaw objects having some state (see
CustomGlenIce) and it's hard to determine whether an object is actually
a derived class, or whether it's a derived class which implements the
same flow law but adds some additional methods, other ways of testing is
error-prone.
 */
class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm,const char prefix[], const NCConfigVariable &conf);
  ~IceFlowLawFactory();
  PetscErrorCode setType(const char[]);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(const char[],
			      PetscErrorCode(*)(MPI_Comm,const char[], const NCConfigVariable &,IceFlowLaw **));
  PetscErrorCode create(IceFlowLaw **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm comm;
  char prefix[256],type_name[256];
  PetscFList type_list;
  const NCConfigVariable &config;
};

// This uses the definition of second invariant from Hutter and several others, namely
// \f$ \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline PetscScalar secondInvariant(PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline PetscScalar secondInvariantDu(const PetscScalar Du[])
{ return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2])); }

#endif /* __materials_hh */

