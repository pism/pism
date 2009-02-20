// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

//! Class containing physical constants and the constitutive relation describing ice.  (Mostly virtual base class.)
class IceType {
public:
  static PetscScalar rho;
  static PetscScalar beta_CC_grad;
  static PetscScalar k;
  static PetscScalar c_p;
  static PetscScalar latentHeat;
  static PetscScalar meltingTemp;
  static PetscScalar n;

  virtual ~IceType() {}
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_FALSE; }
  // this one returns nu; applies to ice shelf/stream approximation
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  // this one returns nu * H; it is adapted to a staggered grid so T1,T2 get averaged
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscInt kbelowH,
                           const PetscInt nlevels, PetscScalar *zlevels,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
  virtual PetscInt integratedStoreSize() const;
  virtual void integratedStore(PetscScalar H, PetscInt kbelowH, PetscInt nlevels, const PetscScalar *zlevels,
                               const PetscScalar T[], PetscScalar store[]) const;
  virtual void integratedViscosity(PetscReal regularization, const PetscScalar store[],
                                   const PetscScalar Du[], PetscScalar *eta, PetscScalar *deta) const;
  virtual PetscScalar exponent() const;
  virtual PetscScalar softnessParameter(const PetscScalar T) const;
  virtual PetscScalar hardnessParameter(const PetscScalar T) const;
};


//! Derived class of IceType for Paterson-Budd (1982)-Glen ice.
class ThermoGlenIce : public IceType {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscInt kbelowH,
                           const PetscInt nlevels, PetscScalar *zlevels,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
  virtual void integratedStore(PetscScalar H, PetscInt kbelowH, PetscInt nlevels, const PetscScalar *zlevels,
                               const PetscScalar T[], PetscScalar store[]) const;
  virtual void integratedViscosity(PetscReal regularization, const PetscScalar store[],
                                   const PetscScalar Du[], PetscScalar *eta, PetscScalar *deta) const;
  virtual PetscScalar softnessParameter(const PetscScalar T) const;
  virtual PetscScalar hardnessParameter(const PetscScalar T) const;
protected:
  static PetscScalar A_cold;  // these four constants from Paterson & Budd (1982)
  static PetscScalar A_warm;
  static PetscScalar Q_cold;
  static PetscScalar Q_warm;
  static PetscScalar crit_temp;
};


//! Derived class of IceType for Hooke (1981)-Glen ice.
class ThermoGlenIceHooke : public ThermoGlenIce {
public:
  virtual PetscScalar softnessParameter(PetscScalar T) const;
protected:
  static PetscScalar A_Hooke;  // these six constants from Hooke (1981)
  static PetscScalar Q_Hooke;
  static PetscScalar C_Hooke;
  static PetscScalar K_Hooke;
  static PetscScalar Tr_Hooke;
  static PetscScalar R_Hooke;
};


//! Derived class of IceType for Arrhenius-Glen ice; \e cold case of Paterson-Budd (1982) ice.
class ThermoGlenArrIce : public ThermoGlenIce {
public:
  virtual PetscScalar softnessParameter(PetscScalar T) const;
  using ThermoGlenIce::flow;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar A() const;  // returns A_cold for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_cold for Paterson-Budd
};


//! Derived class of IceType for Arrhenius-Glen ice; \e warm case of Paterson-Budd (1982) ice.
class ThermoGlenArrIceWarm : public ThermoGlenArrIce {
public:
  virtual PetscScalar A() const;  // returns A_warm for Paterson-Budd
  virtual PetscScalar Q() const;  // returns Q_warm for Paterson-Budd
};



struct GKparts {
  PetscScalar eps_total, eps_diff, eps_disl, eps_basal, eps_gbs;
};

//! Derived class of IceType for a hybrid of Goldsby-Kohlstedt (2001) ice in SIA, with Paterson-Budd (1982)-Glen behavior when needed in viscosity form.
class HybridIce : public ThermoGlenIce {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_TRUE; }
  GKparts flowParts(const PetscScalar stress, const PetscScalar temp,
                    const PetscScalar pressure) const;

protected:
  static PetscScalar  V_act_vol;
  static PetscScalar  d_grain_size;
  //--- diffusional flow ---
  static PetscScalar
  diff_crit_temp, diff_V_m, diff_D_0v, diff_Q_v, diff_D_0b, diff_Q_b, diff_delta;
  //--- dislocation creep ---
  static PetscScalar
  disl_crit_temp, disl_A_cold, disl_A_warm, disl_n, disl_Q_cold, disl_Q_warm;
  //--- easy slip (basal) ---
  static PetscScalar
  basal_A, basal_n, basal_Q;
  //--- grain boundary sliding ---
  static PetscScalar
  gbs_crit_temp, gbs_A_cold, gbs_A_warm, gbs_n, gbs_Q_cold, p_grain_sz_exp, gbs_Q_warm;
};


// Derived class of HybridIce; for testing only.
class HybridIceStripped : public HybridIce {
public:
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure, const PetscScalar gs) const;
  virtual PetscTruth usesGrainSize() const { return PETSC_FALSE; }
protected:
  static PetscScalar d_grain_size_stripped;
};


//! Physical constants describing lithosphere thermal properties.
class BedrockThermalType {
public:
  static PetscScalar rho;
  static PetscScalar k;
  static PetscScalar c_p;
};


//! Physical constants describing lithosphere mechanical properties.
class DeformableEarthType {
public:
  static PetscScalar rho;
  static PetscScalar D;
  static PetscScalar eta;
};


//! Physical constants describing ocean water properties.
class SeaWaterType {
public:
  static PetscScalar rho;
  static PetscScalar beta_CC_grad;
  static PetscScalar defaultShelfBaseMassRate;
};


//! Physical constants describing pure water properties.
class FreshWaterType {
public:
  static PetscScalar rho;
  static PetscScalar beta_CC_grad;
};


// Class containing constitutive relation for till in SIA sliding law; NOT RECOMMENDED.
class BasalTypeSIA {
public:
  virtual PetscScalar velocity(PetscScalar sliding_coefficient,
                               PetscScalar stress);
  virtual ~BasalTypeSIA() {}
};


//! Class containing physical constants and the constitutive relation describing till for SSA.
/*!
This \e pseudo -plastic type can actually describe anything from linearly 
viscous till to purely plastic till.
 */
class PlasticBasalType {
public:
  PlasticBasalType(const PetscScalar regularizationConstant, const PetscTruth pseudoPlastic,
                   const PetscScalar pseudoExponent, const PetscScalar pseudoUThreshold);
  virtual PetscErrorCode printInfo(const int verbthresh, MPI_Comm com);
  virtual PetscScalar drag(const PetscScalar tauc,
                           const PetscScalar vx, const PetscScalar vy);
  // Like drag() but also get the derivative of drag with respect to \f$ alpha=\frac 1 2 \abs{u}^2 \f$.
  virtual void dragWithDerivative(PetscReal tauc, PetscScalar vx, PetscScalar vy, PetscScalar *drag, PetscScalar *ddrag) const;
  virtual ~PlasticBasalType() {} // class w virtual methods needs virtual destructor?

  PetscScalar plastic_regularize, pseudo_q, pseudo_u_threshold;
  PetscTruth  pseudo_plastic;
};


#endif /* __materials_hh */

