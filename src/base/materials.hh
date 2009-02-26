// Copyright (C) 2004-2009 Jed Brown and Ed Bueler
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

/*******************
REGARDING IceType and HybridIce:  The main hierarchy is:
  IceType <- ThermoGlenIce <- HybridIce <- HybridIceStripped,
where "<-" means "derived class".  IceType is a virtual
class; it should never be used "as is".

ThermoGlenIce     means *Paterson-Budd* version of Arhennius relation
ThermoGlenArrIceHooke  means *Hooke* version ...
ThermoGlenArrIce       uses only cold part of Paterson-Budd
ThermoGlenArrIceWarm   uses only warm part of Paterson-Budd
HybridIce         means *Goldsby-Kohlstedt* flow law where vMask=SHEET,
                  and otherwise Paterson-Budd
HybridIceStripped means, where SHEET, G-K without the pressure dependence and without the
                  diffusional part; also grain size fixed at 3mm

Note each IceType has both a forward flow law ("flow") and an
inverted-and-vertically-integrated flow law ("effectiveViscosityColumn").  Only the
former form of the flow law is known for Goldsby-Kohlstedt.  If one can
invert-and-vertically-integrate the G-K law then one can build a "trueGKIce"
derived class.
*******************/

//! Abstract class containing physical constants and the constitutive relation describing ice.
/*!
This is the interface which most of PISM uses for rheology.
 */
class IceType {
public:
  // ideally these would be protected or private:
  PetscScalar rho,          //!< density
              beta_CC_grad, //!< Clausius-Clapeyron gradient
              k,            //!< thermal conductivity
              c_p,          //!< specific heat capacity
              latentHeat,   //!< latent heat capacity
              meltingTemp;  //!< melting temperature

  IceType(MPI_Comm c,const char pre[]);
  virtual ~IceType() {}
  virtual PetscErrorCode setFromOptions() {return 0;}
  virtual PetscErrorCode printInfo(PetscInt) const {return 0;}
  virtual PetscErrorCode view(PetscViewer) const {return 0;}
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const = 0;
  // returns nu * H; it is adapted to a staggered grid so T1,T2 get averaged
  virtual PetscScalar effectiveViscosityColumn(PetscScalar H, PetscInt kbelowH, const PetscScalar *zlevels,
                                               PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y,
                                               const PetscScalar *T1, const PetscScalar *T2) const = 0;
  virtual PetscInt integratedStoreSize() const = 0;
  virtual void integratedStore(PetscScalar H, PetscInt kbelowH, const PetscScalar zlevels[],
                               const PetscScalar T[], PetscScalar store[]) const = 0;
  virtual void integratedViscosity(const PetscScalar store[],const PetscScalar Du[],PetscScalar *nuH, PetscScalar *dnuH) const = 0;
  // This is not a natural part of IceType since it doesn't make any sense for plenty
  // of rheologies.  Nonetheless we need some exponent to compute the coordinate
  // transformation in IceModel::computeDrivingStress (see iMgeometry.cc).
  virtual PetscScalar exponent() const = 0;
  // This is also not a natural part of IceType, but it is needed to invert the Sigma to obtain strain rate in
  // IceModel::correctSigma().  This method can reside here until a plan for generalizing correctSigma has been agreed
  // upon.
  virtual PetscScalar hardnessParameter(PetscScalar T) const = 0;
protected:
  MPI_Comm comm;
  char prefix[256];
};


PetscTruth IceTypeIsPatersonBuddCold(IceType *);
PetscTruth IceTypeUsesGrainSize(IceType *);


class CustomGlenIce : public IceType {
public:
  CustomGlenIce(MPI_Comm c,const char pre[]);
  PetscErrorCode setDensity(PetscReal);
  PetscErrorCode setExponent(PetscReal);
  PetscErrorCode setHardness(PetscReal);
  PetscErrorCode setSoftness(PetscReal);
  PetscErrorCode setSchoofRegularization(PetscReal vel,PetscReal len);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar flow(PetscScalar,PetscScalar,PetscScalar,PetscScalar) const;
  virtual PetscScalar effectiveViscosityColumn(PetscScalar,PetscInt,const PetscScalar[],
                                               PetscScalar,PetscScalar,PetscScalar,PetscScalar,
                                               const PetscScalar[],const PetscScalar[]) const;
  virtual PetscInt integratedStoreSize() const;
  virtual void integratedStore(PetscScalar H, PetscInt kbelowH, const PetscScalar *zlevels,
                               const PetscScalar T[], PetscScalar store[]) const;
  virtual void integratedViscosity(const PetscScalar store[], const PetscScalar Du[], PetscScalar *eta, PetscScalar *deta) const;
  virtual PetscScalar exponent() const;
  virtual PetscScalar hardnessParameter(PetscScalar T) const;
private:
  PetscReal exponent_n,softness_A,hardness_B,schoofVel,schoofLen,schoofReg;
};


//! Derived class of IceType for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceType {
public:
  ThermoGlenIce(MPI_Comm c,const char pre[]);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode printInfo(PetscInt) const;
  virtual PetscErrorCode view(PetscViewer) const;
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;
  virtual PetscScalar effectiveViscosityColumn(PetscScalar,PetscInt,const PetscScalar[],
                                               PetscScalar,PetscScalar,PetscScalar,PetscScalar,
                                               const PetscScalar[],const PetscScalar[]) const;
  virtual PetscInt integratedStoreSize() const;
  virtual void integratedStore(PetscScalar H, PetscInt kbelowH, const PetscScalar *zlevels,
                               const PetscScalar T[], PetscScalar store[]) const;
  virtual void integratedViscosity(const PetscScalar store[], const PetscScalar Du[], PetscScalar *eta, PetscScalar *deta) const;
  virtual PetscScalar exponent() const;
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  virtual PetscScalar hardnessParameter(PetscScalar T) const;
protected:
  PetscReal schoofLen,schoofVel,schoofReg,
            A_cold, A_warm, Q_cold, Q_warm,  // these four constants from Paterson & Budd (1982)
            crit_temp, n;
};


//! Derived class of IceType for Hooke (1981)-Glen ice.  Only changes A(T) factor from ThermoGlenIce.
class ThermoGlenIceHooke : public ThermoGlenIce {
public:
  ThermoGlenIceHooke(MPI_Comm c,const char pre[]);
  virtual PetscScalar softnessParameter(PetscScalar T) const;
protected:
  PetscReal A_Hooke, Q_Hooke, C_Hooke, K_Hooke, Tr_Hooke, R_Hooke; // constants from Hooke (1981)
};


//! Derived class of IceType for Arrhenius-Glen ice; \e cold case of Paterson-Budd (1982) ice.
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  ThermoGlenArrIce(MPI_Comm c,const char pre[]) : ThermoGlenIce(c,pre) {}
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  using ThermoGlenIce::flow;
  virtual PetscScalar flow(PetscScalar,PetscScalar,PetscScalar,PetscScalar) const;
  virtual PetscScalar A() const;  // returns A_cold for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_cold for Paterson-Budd
};


//! Derived class of IceType for Arrhenius-Glen ice; \e warm case of Paterson-Budd (1982) ice.
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  ThermoGlenArrIceWarm(MPI_Comm c,const char pre[]) : ThermoGlenArrIce(c,pre) {}
  virtual PetscScalar A() const;  // returns A_warm for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_warm for Paterson-Budd
};



struct GKparts {
  PetscScalar eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};


//! Derived class of IceType for a hybrid of Goldsby-Kohlstedt (2001) ice in SIA, with Paterson-Budd (1982)-Glen behavior when needed in viscosity form (e.g. SSA).
class HybridIce : public ThermoGlenIce {
public:
  HybridIce(MPI_Comm c,const char pre[]);
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
  HybridIceStripped(MPI_Comm c,const char pre[]);
  virtual PetscScalar flow(PetscScalar stress, PetscScalar temp, PetscScalar pressure, PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_FALSE; }
protected:
  PetscReal d_grain_size_stripped;
};


//! Physical constants describing lithosphere thermal properties.
class BedrockThermalType {
public:
  BedrockThermalType();
  // these should be protected or private
  PetscReal rho,   //!< density
            k,     //!< thermal conductivity
            c_p;   //!< specific heat capacity
};


//! Physical constants describing lithosphere mechanical properties.
class DeformableEarthType {
public:
  DeformableEarthType();
  // these should be protected or private
  PetscReal rho,   //!< density
            D,     //!< lithosphere flexural rigidity
            eta;   //!< half-space (mantle) viscosity
};


//! Physical constants describing ocean water properties.
class SeaWaterType {
public:
  SeaWaterType();
  PetscReal rho,          //!< density
            beta_CC_grad; //!< Clausius-Clapeyron gradient (melting temperature point per meter)
};


//! Physical constants describing pure water properties.
class FreshWaterType {
public:
  FreshWaterType();
  PetscReal rho,          //!< density
            beta_CC_grad; //!< Clausius-Clapeyron gradient (melting temperature point per meter)
};


//! Class containing constitutive relation for till in SIA sliding law; NOT RECOMMENDED.
class BasalTypeSIA {
public:
  virtual PetscScalar velocity(PetscScalar sliding_coefficient,
                               PetscScalar stress);
  virtual ~BasalTypeSIA() {}  // need virtual destructor because has virtual methods
};


//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class PlasticBasalType {
public:
  PlasticBasalType(PetscScalar regularizationConstant, PetscTruth pseudoPlastic,
                   PetscScalar pseudoExponent, PetscScalar pseudoUThreshold);
  virtual PetscErrorCode printInfo(int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(PetscScalar tauc,
                           PetscScalar vx, PetscScalar vy);
  // Also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy, PetscScalar *drag, PetscScalar *ddrag) const;
  virtual ~PlasticBasalType() {} // class w virtual methods needs virtual destructor?

  PetscReal   plastic_regularize, pseudo_q, pseudo_u_threshold;
  PetscTruth  pseudo_plastic;
};


#define ICE_CUSTOM  "custom"        /* Plain isothermal Glen with customizable parameters */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_HYBRID  "hybrid"        /* Goldsby-Kohlstedt for SIA, PB for SSA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

/*!
Excerpt from Jed email to ELB, 2/21:

IceModel creates an IceFactory in it's constructor.  Derived classes
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

To test whether an IceType object is true to a particular constitutive
relation, it's best to just create a reference object and test at a few
values.  Since I envision IceType objects having some state (see
CustomGlenIce) and it's hard to determine whether an object is actually
a derived class, or whether it's a derived class which implements the
same flow law but adds some additional methods, other ways of testing is
error-prone.
 */
class IceFactory {
public:
  IceFactory(MPI_Comm,const char prefix[]);
  ~IceFactory();
  PetscErrorCode setType(const char[]);
  PetscErrorCode setTypeByNumber(int);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(const char[],PetscErrorCode(*)(MPI_Comm,const char[],IceType **));
  PetscErrorCode create(IceType **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm comm;
  char prefix[256],type_name[256];
  PetscFList type_list;
};

// This uses the definition of second invariant from Hutter and several others, namely
// \f$ \frac 1 2 D_{ij} D_{ij} \f$ where incompressibility is used to compute \f$ D_{zz} \f$
static inline PetscScalar secondInvariant(PetscScalar u_x, PetscScalar u_y, PetscScalar v_x, PetscScalar v_y)
{ return 0.5 * (PetscSqr(u_x) + PetscSqr(v_y) + PetscSqr(u_x + v_y) + 0.5*PetscSqr(u_y + v_x)); }

// The second invariant of a symmetric strain rate tensor in compressed form [u_x, v_y, 0.5(u_y+v_x)]
static inline PetscScalar secondInvariantDu(const PetscScalar Du[])
{ return 0.5 * (PetscSqr(Du[0]) + PetscSqr(Du[1]) + PetscSqr(Du[0]+Du[1]) + 2*PetscSqr(Du[2])); }

#endif /* __materials_hh */

